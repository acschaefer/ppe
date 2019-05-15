function runexp
% RUNEXP Perform plane extraction using different algorithms.

%% Set parameters.
% Set maximum squared error increments.
e = (0.02:0.02:1).^2;

% Downsampling for faster test runs.
downsmpl = false;

% Define the algorithms to use for plane extraction.
algorithmname = {'ppe','msac'};

% Create the output directory.
outdir = 'result';
[errorcode,msg] = mkdir(outdir);
if errorcode < 1
    error(['Failed to create output directory ''', outdir, ''': ', msg])
end

% Read datasets.
dataset = load(fullfile('data', 'dataset_test.mat'));
pc = dataset.pc;

% Read PPE parameters.
paramppe = load(fullfile(outdir,'paramppe.mat'));
lmax = paramppe.lmaxopt;

% Read MSAC parameters.
parammsac = load(fullfile(outdir,'parammsac.mat'));
msacdist = parammsac.distopt;
msacperc = parammsac.percopt;

% Downsampling for faster test runs.
if downsmpl
    for i = 1:numel(pc) %#ok<UNRCH>
        ix = 220:230; iy = ix;
        pci = pointCloud(pc{i}.Location(ix, iy, :));
        pci.Intensity = pc{i}.Intensity(ix, iy);
        pc{i} = pci;
    end
end

% Loop over all scans.
disp('Extracting planes from datasets ...')
numpln = numel(e) + 1;
planeres = struct('plane', [], 'steps', [], 'emax', []);
for ids = 1 : size(pc,1)
    % Create result file.
    resultfile = fullfile(outdir,['extrpln_ids',num2str(ids),'.mat']);
    fprintf(['Saving to ', resultfile, '\n']);
    save(resultfile, 'e', 'lmax', 'msacdist', 'msacperc', 'algorithmname');
    
    % Create result structure
    pln = repmat(planeres,[size(pc,2),size(pc,1),numpln,numel(algorithmname)]);
    t = NaN(size(pln));
    try
        % Loop over all pointclouds.
        parfor ipc = 1 : size(pc,2)
            % Loop over all datasets.
            plni = repmat(planeres,[size(pc,1),numpln,numel(algorithmname)]);
            ti = NaN(size(plni));
            % Set lmax
            lmaxi = lmax(min([ids, 2])); %#ok<PFBNS>
            % Add number of ground-truth planes to stopping criteria
            n = max([sum((unique(pc{ids,ipc}.Intensity) >= 10)), 1]);
            % Extract planes using PPE.
            tStart = tic;
            plnres = pcextrpln(pc{ids,ipc}, ...
                'e', e, 'n', n, 'lmax', lmaxi, ...
                'display', 'none', 'device', 'gpu');
            plni(ids,1:numel(plnres),1,1) = plnres;
            ti(ids,1,1,1) = toc(tStart);
            % Extract planes using MSAC.
            msdist = msacdist(min([ids,2])); %#ok<PFBNS>
            msperc = msacperc(min([ids,2])); %#ok<PFBNS>
            tStart = tic;
            plni(ids,1,2,1).plane = extrplnmsac(pc{ids,ipc}, ...
                [msdist, msperc]);
            ti(ids,1,2,1) = toc(tStart);

            % Save results of this iteration.
            pln(ipc,:,:,:) = plni;
            t(ipc,:,:,:) = ti;

            % Update progress display.
            fprintf('Progress: %i/%i.\n', ipc, size(pc,2))
        end
    catch me
        % Save results before processing error.
        saveresult
        rethrow(me)
    end
    % Save results.
    saveresult
end

    function saveresult
        % SAVERESULT Save results to file.
        pln = permute(pln, [2,1,4,3]);
        t = permute(t, [2,1,4,3]);
        save(resultfile, 'pln', 't', '-append')        
    end

end
