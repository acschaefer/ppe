function mergeplnres
%MERGEPLNRES Merge plane extration results from two files.

ids = 1;

% Set file names
file = ['extrpln_ids', num2str(ids), '_old.mat'];
file1 = 'extrplnfeng_segcomp.mat';
if ids > 1
    file1 = 'extrplnfeng_synpeb.mat';
end
fileout = ['extrpln_ids', num2str(ids), '.mat'];

% Print file names
fprintf('%s + %s -> %s\n', file, file1, fileout);

% Load data files
plndata = load(fullfile('result', file));
plndata1 = load(fullfile('result', file1));

% Load data from first file
e = plndata.e;
lmax = plndata.lmax;
algorithmname = plndata.algorithmname;
msacdist = plndata.msacdist;
msacperc = plndata.msacperc;
pln = plndata.pln;
t = plndata.t;

% Load data from second file
algorithmname1 = plndata1.algorithmname;
pln1 = plndata1.pln;
t1 = plndata1.t;

% Merge
algorithmname{end+1} = algorithmname1;
for ipc = 1:size(pln, 2)
    plni = struct('plane', pln1{ids,ipc,3}, 'steps', 0, 'emax', 0);
    pln(ids,ipc,3,1) = plni;
end
t(ids,:,3,1) = t1(ids,:,3);

% Save
outfile = fullfile('result',fileout);
save(outfile, 'e', 'lmax', 'algorithmname', 'msacdist', 'msacperc', 'pln', 't');

fprintf('Done!\n');
end

