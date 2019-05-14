function pc = synpebread(file)
% SYNPEBREAD Read SynPEB dataset.
%   PC = SYNPEBREAD(FILE) reads the SynPEB dataset file whose file name is
%   given by the character vector FILE and returns the corresponding point
%   cloud PC. PC is a pointCloud object whose colors and intensities encode
%   the segmentation ground truth: Points belonging to the same cluster
%   have the same color and intensity.
%
%   Example:
%      gendata
%      pc = synpebread('pcd/data/synpeb/train/gt/pc_01.pcd');
%      pcshow(pc)

% Copyright 2018 Daniel Buescher

pc = pcread(file);
segmentation = cast(pc.Intensity, 'uint32');
max_id = max(reshape(segmentation, size(pc.Intensity,1) * size(pc.Intensity,2), 1));
color_map = randomColorMap(max_id + 1);
segmentation_color = cast(ind2rgb(segmentation, color_map) *255, 'uint8');
pc.Color = segmentation_color;
gtang = [];

end

%% Returns a random color map
function [color_map] = randomColorMap(num_colors)
    color_map = zeros(num_colors,3);
    for i=1:num_colors
        color_map(i,:) = [rand(), rand(), rand()];
    end
end
