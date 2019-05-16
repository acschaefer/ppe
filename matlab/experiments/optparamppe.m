function optparamppe
% OPTPARAMPPE Find optimal parameters for probabilistic plane extraction.

%% Create result file.
% Define the algorithms to use for plane extraction.
algorithmname = 'ppe'; %#ok<NASGU>

% Create the output directory.
outdir = 'result';
[errorcode,msg] = mkdir(outdir);
if errorcode < 1
    error(['Failed to create output directory ''', outdir, ''': ', msg])
end

% Create result file.
resultfile = fullfile(outdir,'paramppe.mat');

%% Optimize plane extraction parameters.
% Read datasets.
dataset = load(fullfile('..', 'data','dataset_train.mat'));
pc = dataset.pc;
gtang = dataset.gtang;
nds = size(pc, 1);

% Set initial parameters and parameter bounds.
x0 = [0.0025,0.03];
xlb = [0.001,0.01];
xub = [0.100,0.50];

% Initialize pattern search output
x = zeros(nds,2);
y = zeros(nds,1);

% Find optimal parameters by evaluating the parameter space.
disp('Optimizing parameters ...')
try
    % Loop over datasets
    for ids = 1:nds
        [x(ids,:), y(ids), ~, ~] = patternsearch( ...
            @(x) evalmsac(pc(ids,:),gtang(ids,:),x), ...
            x0, [], [], [], [], xlb, xub);
    end
catch me
    % Save results before processing error.
    saveresult
    rethrow(me)
end

% Save results.
saveresult

    function saveresult
        % SAVERESULT Save results to file.
        save(resultfile, 'algorithmname', 'x', 'y')
    end

end

function y = evalmsac(pc,gtang,x)
% EVALMSAC Evaluate performance of MSAC plane extraction.

% Compare tolerance
tcomp = 0.8;

% Results
ncorrseg = zeros(numel(pc), 1);

% Loop over all point clouds.
pln = cell(size(pc));
parfor i = 1 : numel(pc)
    % Extract planes from point cloud.
    pln{i} = extrplnmsac(pc{i}, x);
    gti = unique(pc{i}.Intensity);
    gti = gti(gti >= 10);
    [~,ncorrseg(i)] = segcompeval(pc{i}, pln{i}, gtang{i}, tcomp);
end

y = -sum(ncorrseg);
fprintf('Negative quality of MSAC with parameters (%f, %f): %f\n', ...
    x(1), x(2), y);

end
