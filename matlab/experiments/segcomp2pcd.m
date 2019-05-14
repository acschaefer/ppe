% SEGCOMP2PCD Convert SegComp point clouds to PCD files.

% Load datasets.
data = [load(fullfile('pcd','data','planeextract','dataset_train.mat')),...
    load(fullfile('pcd','data','planeextract','dataset_test.mat'))];

% Write training point clouds.
dataname = {'train', 'test'};
for ids = 1 : numel(dataname)
    pc = data(ids).pc(1,:);
    for ipc = 1 : numel(pc)
        % Convert point cloud to struct.
        s = struct('x', pc{ipc}.Location(:,:,1), ...
            'y', pc{ipc}.Location(:,:,2), ...
            'z', pc{ipc}.Location(:,:,3), ...
            'i', pc{ipc}.Intensity);

        % Save point cloud.
        pcdwrite(s, sprintf('segcomp_%s_%02i.pcd',dataname{ids},ipc));
    end
end
