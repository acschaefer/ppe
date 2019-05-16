% PLOTPROB Create plots to demonstrate problems with SegComp dataset.

% Create data.
pc = segcompread(fullfile('..', 'data','segcomp','perc.test.23'));
[c,~,i] = unique(pc.Intensity(:));
pc.Color = ind2rgb8(reshape(i,size(pc.Location(:,:,1))),lines);
seg = permute(pc.Color,[2,1,3]);

% Display data.
figure('InvertHardcopy', 'off')
pcshow(pc, 'MarkerSize', 30);
campos([-1.4266, 0.8919, -0.6199])
figure
imshow(seg)
imwrite(reshape(i,size(pc.Intensity))', lines(numel(c)), 'plotprob.png')
