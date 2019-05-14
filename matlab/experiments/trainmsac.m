function trainmsac
% TRAINMSAC Find optimal parameters for MSAC plane extraction.

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

% Create grid of parameters.
dist = [0.0005; 0.002] .* (1:40);
perc = 0.01 * (1:30);

% Initialize search output.
nds = numel(dsidx);
ndist = size(dist, 2);
y = NaN([nds, ndist, numel(perc)]);
sy = size(y);
n = numel(y);

% Evaluate all combinations of parameters.
fprintf('Number of jobs: %i\n', n);
parfor i = 1 : n
    [ids,idist,iperc] = ind2sub(sy, i);
    fprintf('Starting i (ids idist iperc): %i (%i %i %i)\n', ...
        i, ids, idist, iperc);
    y(i) = evalmsac(data.pc(dsidx(ids),:), ...
        [dist(ids,idist),perc(iperc)]);
end

% Determine the optimal combination of parameters.
distopt = NaN(nds, 1);
percopt = NaN(nds, 1);
for ids = 1 : nds
    % Determine the indices of the optimal parameters.
    yds = y(ids,:,:);
    [ymin,imin] = min(yds(:));
    [idist,iperc] = ind2sub([ndist,numel(perc)], imin);
    
    % Determine the values of the optimal parameters.
    distopt(ids) = dist(ids,idist);
    percopt(ids) = perc(iperc);
    
    % Show the optimal parameters on the console.
    fprintf('Dataset %s:\n', data.datasetname{dsidx(ids)})
    fprintf('\tMaximum fraction of detected planes: %f.\n', -ymin);
    fprintf('\tCorresponding parameters: dist=%.3f, perc=%.2f.\n', ...
        distopt(ids), percopt(ids));
    
end

%% Save results.
algorithmname = 'MSAC'; %#ok<NASGU>
save(fullfile(outdir,'parammsac.mat'), ...
    'algorithmname', 'dist', 'perc', 'y', 'distopt', 'percopt')

%% Plot results.
try
    for ids = 1 : nds
        subplot(nds,1,ids);
        surf(perc, dist(ids,:), squeeze(-y(ids,:,:)));
        xlabel('percentage');
        ylabel('distance');
        zlabel('fraction of correctly segmented planes')
        title(data.datasetname{dsidx(ids)});
    end
catch
    warning('Plotting failed.');
end

end

function y = evalmsac(pc, x)
% EVALMSAC Evaluate performance of MSAC plane extraction.

% Define required overlap for a plane to be classified as correctly
% detected.
tcomp = 0.8;

% Preallocate results.
cs = NaN(size(pc));

% Loop over all point clouds.
warning('off', 'vision:ransac:maxTrialsReached')
for i = 1 : numel(pc)
    
    % Downsampling for faster test runs.
    if false
        ix = 220:260; iy = ix;
        pci = pointCloud(pc{i}.Location(ix, iy, :));
        pci.Intensity = pc{i}.Intensity(ix, iy);
        pc{i} = pci;
    end
    
    % Extract planes from point cloud.
    pln = extrplnmsac(pc{i}, x);
    
    % Compute the fraction of correctly segmented planes.
    [~,ncorrseg] = segcompeval(pc{i}, pln, [], tcomp);
    cs(i) = ncorrseg / numel(unique(pc{i}.Intensity));
end

% Use mean fraction of correctly segmented planes as performance metric.
y = -mean(cs);

end
