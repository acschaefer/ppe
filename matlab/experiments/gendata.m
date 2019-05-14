function gendata
% GENDATA Download dataset to evaluate plane extraction.

%% Create output directory.
% Define directory where to write the dataset files.
outdir = fullfile('pcd','data','planeextract');

% Create the output directory.
[errorcode,msg] = mkdir(outdir);
if errorcode < 1
    error(['Failed to create output directory ''', outdir, ''': ', msg])
end

%% Read SegComp dataset.
% Define the URL from where to download the SegComp dataset.
url = 'figment.csee.usf.edu';

% Define the path to the SegComp dataset on the server.
srvfolder = 'pub/segmentation-comparison/';

% Specify the path to the dataset on the local machine.
datafolder = fullfile(fileparts(mfilename('fullpath')), '..', '..', ...
    'data', 'segcomp');
tarfolder = fullfile(datafolder, srvfolder);
zfolder = fullfile(datafolder, 'z');

% Download test and training part of the dataset.
pctrain = {};
gtangtrain = {};
pctest = {};
gtangtest = {};
dstype = {'test','train'};
for dst = dstype
    % Define the name of the dataset file.
    file = ['PERC-', upper(dst{:}), '-IMAGES.tar'];

    % Download the dataset only if not already downloaded.
    if ~exist(fullfile(tarfolder,file), 'file')
        % If the folder does not exist, create it.
        if ~exist(tarfolder, 'dir')
            mkdir(tarfolder);
        end

        % Connect to the FTP server.
        if ~exist('ftpsrv', 'var')
            ftpsrv = ftp(url);
            cleaner = onCleanup(@() ftpsrv.close);
        end
        
        % Download the compressed dataset file.
        mget(ftpsrv, [srvfolder,file], datafolder);
    end

    % Untar the dataset to a set of Z files.
    zfile = untar(fullfile(tarfolder,file), zfolder);

    % Uncompress the Z files.
    for i = find(endsWith(zfile, '.z', 'IgnoreCase', true))
        % Define system command to uncompress the file.
        if ispc % Windows.
            % Use the version of 7-zip shipped with MATLAB to uncompress
            % the file.
            cmd = ['7z.exe x "',zfile{i},'" -y -o"',datafolder,'" > NUL'];
        elseif isunix % Linux or Mac.
            [~,name,ext] = fileparts(zfile{i});
            cmd = ['cd ' , datafolder, ' && cp ', zfile{i}, ...
                ' . && uncompress -f ', [name,ext]];
        else
            error(['Unsupported platform. ', upper(mfilename), ...
                ' only supports Windows, Unix, and Mac.'])
        end

        % Uncompress the file.
        if system(cmd) ~= 0
            error('Failed to extract Z files.')
        end
    end

    % Determine the names of the extracted data files.
    datafile = {};
    for file = zfile
        if regexp(file{:}, ['perc.', dst{:}, '.\d+[.]Z$'])
            [~,name] = fileparts(file{:});
            datafile{1,end+1} = fullfile(datafolder, name); %#ok<AGROW>
        end
    end

    % Read the SegComp dataset.
    for file = datafile
        switch dst{:}
            case 'train'
                [pctrain{1,end+1},gtangtrain{1,end+1}] ...
                    = segcompread(file{:}, 'yrange', [1,508]); %#ok<AGROW>
            case 'test'
                [pctest{1,end+1},gtangtest{1,end+1}] ...
                    = segcompread(file{:}, 'yrange', [1,508]); %#ok<AGROW>
        end
    end
end

%% Read SynPEB dataset.
% Define the URL from where to download the SynPEB dataset.
url = 'http://synpeb.cs.uni-freiburg.de/synpeb.zip';

% Specify the path to the dataset on the local machine.
datafolder = fullfile(fileparts(mfilename('fullpath')), '..', '..', ...
    'data', 'synpeb');
    
% Download and uncompress the dataset.
unzip(url, datafolder);

% Loop over training and test dataset.
for dst = dstype
    % Loop over all folders with different variances.
    ids = 2;
    for varfolder = dir(fullfile(datafolder, dst{:}))'
        % Loop over all PCD files in the folder.
        ipc = 1;
        for file = dir(fullfile(varfolder.folder, varfolder.name))'
            % Make sure the file is a PCD file.
            if endsWith(fullfile(file.folder,file.name), '.pcd', ...
                    'IgnoreCase', true)
                % Add the point cloud to the dataset.
                f = fullfile(file.folder, file.name);
                switch dst{:}
                    case 'train'
                        pctrain{ids,ipc} = synpebread(f); %#ok<AGROW>
                        gtangtrain{ids,ipc} = []; %#ok<AGROW>
                    case 'test'
                        pctest{ids,ipc} = synpebread(f); %#ok<AGROW>
                        gtangtest{ids,ipc} = []; %#ok<AGROW>
                end
                ipc = ipc + 1;
            end
        end
        ids = ids + (ipc>1);
    end
end
    
%% Save datasets.
% Define the names of the datasets.
datasetname = {'SegComp'; 'SynPEB 0.0'; 'SynPEB 0.5'; 'SynPEB 1.0'; ...
    'SynPEB 2.0'; 'SynPEB 4.0'}; %#ok<NASGU>

% Write the point clouds to MAT files.
pc = pctrain; %#ok<NASGU>
gtang = gtangtrain; %#ok<NASGU>
save(fullfile(outdir, 'dataset_train.mat'), 'pc', 'gtang', 'datasetname');
pc = pctest; %#ok<NASGU>
gtang = gtangtest; %#ok<NASGU>
save(fullfile(outdir, 'dataset_test.mat'), 'pc', 'gtang', 'datasetname');

end
