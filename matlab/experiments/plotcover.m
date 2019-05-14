function plotcover
% PLOTCOVER Create cover image for plane extraction paper.

% Read point cloud.
pc = pcread(fullfile('pcd','data','planeextract','cover.pcd'));

% Extract planes.
res = pcextrpln(pc, 'e', (0.24-0.1:0.02:0.24+0.1).^2, 'lmax', 0.18, ...
    'device', 'gpu');

% Save result.
save(fullfile('pcd','result','planeextract','cover.mat'), 'pc', 'res');

% Display result.
plotseg(pc, res(6).plane)

end
