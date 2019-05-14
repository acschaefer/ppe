function trainppe
% TRAINPPE Find optimal parameters for probabilistic plane extraction.

%% Create output directory.
outdir = fullfile('pcd','result','planeextract');
[errorcode,msg] = mkdir(outdir);
if errorcode < 1
    error(['Failed to create output directory ''', outdir, ''': ', msg])
end

%% Optimize parameters.
% Read training datasets.
data = load(fullfile('pcd','data','planeextract','dataset_train.mat'));
dsidx = [1,5];

% Create parameter grid.
lmax = [0.004:0.004:0.04; 0.02:0.02:0.2];
e = (0.02:0.02:1).^2;

% Initialize search output.
nlmax = size(lmax, 2);
nds = numel(dsidx);
npc = size(data.pc, 2);
ne = numel(e);
y = NaN([nds*nlmax*npc, ne]);

% Evaluate all combinations of parameters.
nit = nds*nlmax*npc;
fprintf('Number of jobs: %i\n', nit);
parfor it = 1 : nit
    [ids, ilmax, ipc] = ind2sub([nds, nlmax, npc], it);
    %if ids ~= 1; continue; end
    fprintf('Starting it (ids ilmax ipc): %i (%i %i %i)\n', ...
        it, ids, ilmax, ipc);
    y(it,:) = evalppe(data.pc{dsidx(ids),ipc}, ...
        e, lmax(ids,ilmax));
end

% Recover proper dimensions
y = reshape(y, [nds, nlmax, npc, ne]);

% Use mean fraction of correctly segmented planes as performance metric.
y = mean(y, 3);

% Determine the optimal combination of parameters.
lmaxopt = NaN(nds, 1);
eopt = NaN(nds, 1);
for ids = 1 : nds
    % Determine the indices of the optimal parameters.
    yds = y(ids,:,:);
    [ymin,imin] = min(yds(:));
    [ilmax,ie] = ind2sub([nlmax,ne], imin);
    
    % Determine the values of the optimal parameters.
    lmaxopt(ids) = lmax(ids,ilmax);
    eopt(ids) = e(ie);
    
    % Show the optimal parameters on the console.
    fprintf('Dataset %s:\n', data.datasetname{dsidx(ids)})
    fprintf('\tMaximum fraction of detected planes: %f.\n', -ymin);
    fprintf('\tCorresponding parameters: lmax=%.3f, e=%f.\n', ...
        lmaxopt(ids), eopt(ids));
    
end

%% Save results.
algorithmname = 'PPE'; %#ok<NASGU>
save(fullfile(outdir,'paramppe.mat'), ...
    'algorithmname', 'lmax', 'e', 'y', 'lmaxopt', 'eopt')

%% Plot results.
try
    for ids = 1 : nds
        subplot(nds,1,ids);
        surf(sqrt(e), lmax(ids,:), squeeze(-y(ids,:,:)));
        xlabel('sqrt(e)');
        ylabel('lmax');
        zlabel('fraction of correctly segmented planes')
        title(data.datasetname{dsidx(ids)});
    end
catch
    fprintf('Plotting failed.\n');
end

end

function y = evalppe(pc, e, lmax)
% EVALMSAC Evaluate performance of probabilistic plane extraction.

% Define required overlap for a plane to be classified as correctly
% detected.
tcomp = 0.8;

% Preallocate results.
cs = NaN(1,numel(e));

% Downsampling for faster test runs.
if false
    ix = 220:240; iy = ix;
    pci = pointCloud(pc.Location(ix, iy, :));
    pci.Intensity = pc.Intensity(ix, iy);
    pc = pci;
end

% Extract planes from point cloud.
pln = pcextrpln(pc, 'e', e, 'lmax', lmax, 'device', 'gpu');

% Compute the fraction of correctly segmented planes.
for ie = 1 : numel(e)
    % Determine parameter index for given error value
    j = 1;
    for j1 = 1:numel(pln)
        if pln(j1).emax <= e(ie)
            j = j1;
        end
    end
    plnj = rmatompln(pln(j));
    [~,ncorrseg] = segcompeval(pc, plnj.plane, [], tcomp);
    gti = unique(pc.Intensity);
    gti = gti(gti >= 10);
    ngt = numel(gti);
    cs(ie) = ncorrseg / ngt;
end

y = -cs;

end
