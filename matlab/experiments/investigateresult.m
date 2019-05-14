% Load dataset.
dataset = load(fullfile('pcd','data','planeextract','dataset_test.mat'));
pc = dataset.pc;
gtang = dataset.gtang;

% Load experimental results.
plndata = load(fullfile('pcd','result','planeextract','extrpln.mat'));
algorithmname = plndata.algorithmname;
pln = plndata.pln;


ipc = 1; ipm = 8; igt1 = 10; igt2 = 13; ims1 = 6; ims2 = 7;
%ipc = 1; ipm = 8; igt1 = 12; igt2 = 13; ims1 = 5; ims2 = 7;

pc = pc{1,ipc}
gtang = gtang{1, ipc}
pln = pln(1,ipc,1,ipm)

plnsi = (pc.Intensity==igt1 | pc.Intensity==igt2);

np = numel(pc.Location)/3;
loc = reshape(pc.Location, [np, 3]);
loc = loc(plnsi,:);
pc1 = pointCloud(loc);

figure;
pcshow(pc1);
campos([0 0.2 0]);

figure;
plotplane(pc, pln.plane([ims1, ims2]));
campos([0 0.2 0]);

