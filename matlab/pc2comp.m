% PC2COMP Compress PCD files in current directory.

d = dir('.');
parfor i = find(endsWith({d.name}, '.pcd', 'IgnoreCase', true))
    file = fullfile(d(i).folder, d(i).name);
    pcwrite(pcread(file), file, 'Encoding', 'compressed');
end
