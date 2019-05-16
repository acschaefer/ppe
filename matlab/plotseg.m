function plotseg(pc,pln)
% PLOTSEG Plot segmentation results of function PCEXTRPLN.

% Determine the plane index of each point in the point cloud.
pi = zeros(pc.Count, 1);
for ipln = 1 : numel(pln)
    if numel(pln(ipln).index) >= 4
        pi([pln(ipln).index]) = ipln;
    end
end
pi = reshape(pi, size(pc.Location,1), size(pc.Location,2));

% Create a random colormap.
ngt = sum((unique(pc.Intensity) >= 10));
rng(0);
cmgt = rand(ngt,3);
cmms = rand(numel(pln),3);

[~, ~, ~, ~, ~, ~, ~, ~, ~, corrseg] = segcompeval(pc, pln, [], 0.8);
[i, j] = find(corrseg);
cmgt(i,:) = cmms(j,:);
cmgt(end+1,:) = [0 0 0];

% Create the segmentation image.
seg = zeros(size(pc.Location));
for ix = 1 : size(seg,1)
    for iy = 1 : size(seg,2)
        if pi(ix,iy) > 0
            seg(ix,iy,:) = cmms(pi(ix,iy),:);
        end
    end
end

% Plot the ground truth segmentation.
subplot(1,2,1)
[un,~,ic] = unique(pc.Intensity);
ic = ic - (numel(un) - ngt);
ic(ic <= 0) = ngt + 1;
imgt = fliplr(reshape(cmgt(ic,:), size(pc.Location)));
seggt = imshow(imgt);
axis equal
axis off
title('Ground truth segmentation')

% Plot the estimated segmentation.
subplot(1,2,2)
imest = fliplr(seg);
segest = imshow(imest);
axis equal
axis off
title('Estimated segmentation')
linkaxes([seggt.Parent,segest.Parent], 'xy')

% Save plots to file.
%imwrite(uint8(imgt*255), 'seggt.png')
%imwrite(uint8(imest*255), 'segest.png')

end
