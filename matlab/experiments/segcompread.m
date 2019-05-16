function [pc,gtang] = segcompread(file, varargin)
% SEGCOMPREAD Read SegComp dataset.
%   PC = SEGCOMPREAD(FILE) reads the SegComp dataset file whose file name
%   is given by the character vector FILE and returns the corresponding
%   point cloud PC. PC is a pointCloud object whose colors and intensities
%   encode the segmentation ground truth: Points belonging to the same
%   cluster have the same color and intensity.
%
%   [PC,GTANG] = SEGCOMPREAD(FILE) additionally returns the angles GTANG
%   between all pairs of ground truth planes.
%
%   SEGCOMPREAD can be further configured by calling
%   SEGCOMPREAD(...,Name,Value) using the following name-value pairs:
%      'xrange'  - Range of columns to select in the image,
%                  i.e. [xmin, xmax]. Defaults to [].
%      'yrange'  - Range of rows to select in the image,
%                  i.e. [ymin, ymax]. Defaults to [].
%
%   Example:
%      gendata
%      pc = segcompread(fullfile('..', 'data', 'segcomp', 'perc.test.0'))
%      pcshow(pc)
%
%   See also POINTCLOUD, PCREAD.

% Copyright 2018 Johan Vertens, Alexander Schaefer, Daniel Buescher

%% Validate input.
% Check the file name.
validateattributes(file, {'char'}, {'scalartext'}, '', 'FILENAME')

% Parse name-value pair arguments.
parser = inputParser;
parser.addParameter('xrange', []);
parser.addParameter('yrange', []);
parser.parse(varargin{:});
xrange = parser.Results.xrange;
yrange = parser.Results.yrange;

%% Read point cloud and segmentation map.
% Check if the files containing the point cloud and the ground truth
% segmentation exist.
pcfile = file;
segfile = [file, '.gt-seg'];
for f = {pcfile, segfile}
    assert(exist(f{:},'file')==2, ['File ',f{:},' does not exist.'])
end

%% Read point cloud and segmentation map.
% Read point cloud and segmentation map.
[pc, ind] = binread(pcfile, xrange, yrange);
seg = imread(segfile)';
seg = seg(ind)';

% Create the colormap.
useg = unique(seg);
jetmap = jet(numel(useg));
rng(0)
colormap(useg,:) = jetmap(randperm(numel(useg)),:);

% Create point cloud with segmentation encoded in point colors and
% intensities.
pc.Color = uint8(reshape(colormap(seg(:),:), [size(seg),3]) * 255);
pc.Intensity = single(seg);

% Read angles between pairs of ground truth planes.
if nargout > 1
    gtangfile = [file, '.gt-ang'];
    if exist(gtangfile,'file') == 2  
        gtang = readGTAngles(gtangfile);
    else
        gtang = [];
        warning(['File ', gtangfile, ' does not exist.'])
    end
end

end

function [pc, ind] = binread(file, xrange, yrange)
% BINREAD Read a pointcloud in SegComp binary format.

%% Open file.
% Try to open the file.
fid = fopen(file, 'r');
if fid == -1
    error(['Cannot read ', file, '.'])
end

% Make sure the file is closed upon function termination.
cleaner = onCleanup(@() fclose(fid));

%% Read header.
header = {};
for c = 1:10
    line = fgets(fid);
    header = [header, line]; %#ok<AGROW>
end
% Read extra bytes.
fgets(fid);

%% Read range image.
rows = str2double(cell2mat(strtrim(extractAfter(header(1), ':'))));
columns = str2double(cell2mat(strtrim(extractAfter(header(2), ':'))));
data = fread(fid, rows*columns,  'int32', 0,  'b');
rangeimage = reshape(data, rows, columns);
rangeimage = cast(rangeimage, 'int32');
rangeimage = bitand(rangeimage, 4095, 'int32');
rangeimage = cast(rangeimage, 'double');

% Determine indices for range selection
xmin = 1;
xmax = columns;
ymin = 1;
ymax = rows;
if numel(xrange) == 2
    xmin = max([xrange(1), xmin]);
    xmax = min([xrange(2), xmax]);
end
if numel(yrange) == 2
    ymin = max([yrange(1), ymin]);
    ymax = min([yrange(2), ymax]);
end
ind = (columns)*((ymin:ymax)-1)'+(xmin:xmax);
columns = xmax - xmin + 1;
rows = ymax - ymin + 1;

% Select images ranges
rangeimage = rangeimage(ind);
reshape(ind, rows, columns);

cartmap = sc2cart(rangeimage);
% Scale unit to [m].
cartmap = cartmap ./ 1e3;
% Swap x and y axis in image
cartmap = permute(cartmap, [2 1 3]);
pc = pointCloud(cartmap);

end

function gtang = readGTAngles(file)
% READGTANGLES Read ground-truth angles.

%% Open file.
% Try to open the file.
fid = fopen(file, 'r');
if fid == -1
    error(['Cannot read ', file, '.'])
end

% Make sure the file is closed upon function termination.
cleaner = onCleanup(@() fclose(fid));

%% Read file.
% Read header
nlines = str2num(fgetl(fid)); %#ok<ST2NM>
% Read GT plane ids and corresponding angles
gtang = NaN(nlines, 3);
for i = 1:nlines
    gtang(i,:) = str2num(fgetl(fid)); %#ok<ST2NM>
end

end

function [x] = sc2cart(rangeimage)
% SC2CART Convert SegComp depth-map to Cartesian format.

h_1=3.0;			% dist (y) between rotating mirror axis
                    % and the parallel laser beam. 
h_2=5.5; 			% dist (y) between nodding mirror axis
                    % and rotating mirror laser intersection.
gamma=45.0*(pi/180.0);	% slope of rotating mirror
theta=45.0*(pi/180.0);	% slope of nodding mirror in mid position 
alpha_0 = 0.0;
beta_0 = 0.0;
H=51.65;
V=36.73;
r_0=830.3;
delta=0.20236;

nr = size(rangeimage,1);
nc = size(rangeimage,2);
x = zeros(nr, nc, 3);

for r=0:nr-1
    for c=0:nc-1
        alpha = alpha_0 + H*(255.5 - c)/nc;
        alpha = alpha*pi/180.0;
        beta  = beta_0  + V*(255.5 - r)/nr;
        beta  = beta*pi/180.0;

        % (dx,dy,dz): where the laser leaves the nodding mirror */

        dz = -(h_1*(1.0-cos(alpha))/tan(gamma)); % Motion of the laser
                                                 % intersection
                                                 % point with the 
                                                 % rotating mirror 
                                                 % as a function of
                                                 % alpha. */

        dy = dz*tan(theta + 0.5*beta);  % Vertical offset of nodding
                                        % mirror intersection due
                                        % to dz */

        dx = (h_2 + dy)*tan(alpha);     % Horizontal offset of laser
                                        %intersection
                                        %point with nodding mirror due
                                        %to alpha */


        % R  : true range: laser light path length (laser travel) from
        % laser emission point to target.
        % ra : range value returned by camera.
        % r_0: (standoff distance) length of laser ray to point for
        % which r=0
        % r_1: laser travel (z) from laserdiode to intersection with
        % rotating
        % mirror: r1 = d1 + dz
        % r_2: laser travel (vector length) from rotating mirror
        % to nodding mirror.
        % r_3: laser travel (vector length) from nodding mirror to
        % target.
        % Thus: R = r + r0 = r1 + r2 + r3

        ra  = double(rangeimage(r+1,c+1)); % range returned by camera

        r_1 = (dz - h_2)/delta;
        r_2 = sqrt(dx*dx + (h_2+dy)*(h_2+dy))/delta;
        r_3 = (ra + r_0 - (r_1 +r_2))*delta;

        x(r+1,c+1,1) = dx + r_3 * sin(alpha);
        x(r+1,c+1,2) = dy + r_3 * cos(alpha)*sin(beta);
        x(r+1,c+1,3) = dz + r_3 * cos(alpha)*cos(beta); 
    end
end

end
