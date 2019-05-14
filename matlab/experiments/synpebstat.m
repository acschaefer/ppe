function synpebstat
% SYNPEBSTAT Compute statistics of SynPEB dataset.

% Read point clouds.
data = load(fullfile('pcd','data','planeextract','dataset_test.mat'));
pc = data.pc(2,:);

% Compute statistics.
np = NaN(size(pc));
npln = NaN(size(pc));
for i = 1 : numel(pc)
    np(i) = pc{i}.Count;
    npln(i) = sum(unique(pc{i}.Intensity(:))>=10);
end

fprintf(['Number of points per scan:', ...
    repmat(' %i',size(np)), '.\n'], np)
fprintf(['Number of planes per scan:', ...
    repmat(' %i',size(np)), '.\n'], npln)
fprintf('Total number of planes: %i.\n', sum(npln))
