function res = pcextrpln(pc, varargin)
% PCEXTRPLN Extract planes from point cloud.
%   RES = PCEXTRPLN(PC, STOP) extracts planes from the point cloud PC and
%   stores the results in the structure RES. PCEXTRPLN uses a probabilistic
%   greedy search algorithm that creates, expands, and merges planes until
%   a given stopping criterion is met. In each step, the algorithm takes
%   the action that increases the overall error the least, i.e. the action
%   that maximizes the probability that the point cloud is in fact caused
%   by the extracted planes. When computing the error, PCEXTRPLN assumes
%   that the sensor provides perfect angluar information and exhibits
%   normally distributed noise in radial direction only.
%
%   PC is a pointCloud object that represents an organized point cloud. The
%   points are assumed to be specified with respect to the sensor
%   coordinate frame. If the points are specified in a frame other than the
%   sensor coordinate frame, they need to be transformed to the sensor
%   coordinate frame before being passed to PCEXTRPLN.
%   
%   STOP consists of up to two name-value pairs that can define an
%   arbitrary number of stopping criteria. The algorithm stops when the
%   least restrictive stopping criterion is met. Whenever the algorithm
%   passes a more restrictive stopping criterion on the way, the
%   corresponding intermediate result is stored in RES. The stopping
%   criteria are:
%      'n'    - number of planes, specified as integers. In the beginning,
%               each point of the point cloud defines a plane. Plane
%               extraction reduces the number of planes to the given
%               values.
%      'e'    - maximum admissible error increment in a plane extraction
%               step, specified as positive real numbers. The algorithm
%               stops if the error increment corresponding to the next
%               creating, extending, or merging step would exceed this
%               value. The error is defined as the sum of the squared
%               distances between the plane and the points represented by
%               the plane. The distance between a point and the plane is
%               measured along the ray the point is located on. The lower
%               'e', the higher the number of planes created, and the
%               higher the accuracy by which the planes approximate the
%               point cloud. The higher 'e', the lower the number of planes
%               created, and the lower the accuracy.
%   If no stopping criterion is given, 'e' is set to the default of 0.01.
%
%   RES is a struct array that contains the results of the plane extraction
%   process. Each struct contains a snapshot of the plane extraction
%   process at the point a stopping criterion was encountered. RES consists
%   of the following fields:
%      'steps' - number of executed plane extraction steps, specified as
%                integer scalar.
%      'emax'  - maximum admissible error increment, specified as real 
%                scalar.
%      'plane' - struct array of all extracted planes that represent at
%                least four points. Each element of 'plane' contains the
%                following fields:
%                   'index' - vector that contains the linear indices into
%                             the points of PC that are represented by the
%                             plane.
%                   'param' - 3x3 matrix whose rows represent the support
%                             vectors of the plane with respect to the
%                             sensor coordinate system.
%                   'error' - real scalar that contains the error incurred
%                             by approximating the points belonging to the
%                             plane by the plane.
%
%   PCEXTRPLN can be further configured by calling
%   PCEXTRPLN(...,Name,Value) using the following name-value pairs:
%      'lmax'    - real scalar that defines the maximum admissible
%                  Cartesian distance between neighboring ray endpoints
%                  that belong to the same plane.
%                  Defaults to Inf.
%      'display' - defines how to visualize the plane extraction process:
%                     'none'  - no visualization.
%                     'iter'  - visualizes each plane creation, extraction,
%                               and merging step.
%                     'final' - shows the plane extraction result only.
%                  Defaults to 'none'.
%      'device'  - PCEXTRPLN can run both on the CPU and on the GPU. The
%                  device can be specified by
%                     'mat'   - Native MATLAB CPU execution.
%                     'cpu'   - CPU execution using Cpufit library.
%                     'gpu'   - GPU execution using Gpufit library.
%                               GPU execution can reduce computation
%                               times significantly. Requires CUDA
%                               compatible graphics card.
%                  Defaults to 'mat'.
%
%   See also LASERSCAN2.EXTRLIN, LINXPLN.

% Copyright 2018 Alexander Schaefer

%% Validate input.
% Check the data type of the point cloud.
validateattributes(pc, {'pointCloud'}, {'scalar'}, '', 'PC')

% Make sure the point cloud is organized.
if ismatrix(pc.Location)
    error('PC must be organized.')
end

% Make sure the size of the point cloud is at least 2x2.
spc = [size(pc.Location,1), size(pc.Location,2)];
if any(spc < 2)
    error('PC must be at least of size 2x2.')
end

% Parse name-value pair arguments.
parser = inputParser;
parser.addParameter('n', []);
parser.addParameter('e', []);
parser.addParameter('lmax', Inf, @(lmax) ...
    validateattributes(lmax, {'numeric'}, ...
    {'real', 'finite', 'positive', 'scalar'}, '', '''lmax'''))
parser.addParameter('display', 'none', @(display) ...
    ~isempty(validatestring(display, {'none','iter','final'})))
parser.addParameter('device', 'mat', @(visualize) ...
    ~isempty(validatestring(visualize, {'mat','cpu','gpu'})))
parser.parse(varargin{:})
nstop = sort(unique(parser.Results.n),'descend');
estop = unique(parser.Results.e);
lmax = parser.Results.lmax;
display = parser.Results.display;
device = parser.Results.device;

% Check stopping criteria.
if ~isempty(nstop)
    validateattributes(nstop, {'numeric'}, ...
        {'integer','positive','finite'}, '', '''n''')
end
if ~isempty(estop)
    validateattributes(estop, {'numeric'}, ...
        {'real','positive','finite'}, '', '''e''')
end

% Set the default stopping criterion.
if isempty(nstop) && isempty(estop)
    estop = 0.01;
end

% Check if cpufit is available.
if strcmpi(device, 'cpu') || strcmpi(device, 'gpu')
    assert(exist('CpufitMex','file')~=0, ['CpufitMex file not found. ', ...
        'Please build Cpufit and add the resulting files to the path.']);
end

% Check if CUDA is available.
if strcmpi(device, 'gpu')
    assert(exist('GpufitMex','file')~=0, ['GpufitMex file not found. ', ...
        'Please build Gpufit and add the resulting files to the path.']);
    assert(gpufit_cuda_available, 'CUDA is not available.')
end

%% Initialize planes.
% Compute the radii and normalized direction vectors of all rays.
np = pc.Count;
l = double(reshape(pc.Location, np, 3, 1));
r = vecnorm(l, 2, 2);
v = l ./ r;

% Handle no-return rays.
v(r==0,:) = NaN;
r(r==0) = NaN;

% For each point, determine the indices of two neighboring points. The
% corresponding three ray direction vectors must by linearly independent.
% The set of the direction vectors, together with the lengths of the three
% rays, will define the orientation of the plane originating in the point.
[sx,sy] = ndgrid(1:spc(1), 1:spc(2));
s = sub2ind(spc, [sx(:), sx(:)+(-1).^(sx(:)>spc(1)/2), sx(:)], ...
    [sy(:), sy(:), sy(:)+(-1).^(sy(:)>spc(2)/2)]);

% Initialize the set of planes.
plncrt = cell(np,1);
plno = cell(np,1);
plne = NaN(np,1);
plnx = r(s);

%% Initialize plane seeds.
% Define the different initial plane configurations for one point.
c = {[0,0,0; 0,1,1; 0,1,1];
    [0,1,0; 0,1,1; 0,1,0];
    [0,1,0; 1,1,0; 0,1,0];
    [0,0,0; 1,1,1; 0,1,0];
    [0,1,0; 1,1,1; 0,0,0];
    [0,0,1; 0,1,1; 0,1,0];
    [1,0,0; 1,1,0; 0,1,0];
    [0,0,0; 1,1,0; 0,1,1];
    [0,1,1; 1,1,0; 0,0,0];
    [0,1,1; 0,1,0; 0,1,0];
    [1,1,0; 0,1,0; 0,1,0];
    [0,0,0; 1,1,1; 1,0,0];
    [0,0,0; 1,1,1; 0,0,1];
    [1,0,0; 1,1,1; 0,0,0];
    [0,0,1; 1,1,1; 0,0,0];
    [0,1,0; 0,1,0; 0,1,1];
    [0,1,0; 0,1,0; 1,1,0]};

% Determine the index addends corresponding to each configuration in the
% organized point cloud.
nc = numel(c);
[dcx,dcy] = cellfun(@find, c, 'UniformOutput', false);
dcx = reshape(cell2mat(dcx),4,nc)' - 2;
dcy = reshape(cell2mat(dcy),4,nc)' - 2;
        
% Initialize the set of plane seeds.
crte = NaN(np,nc);
crto = cell(np,nc);
crtx = cell(np,nc);

% Determine the Cartesian distances between the points in vertical and
% horizontal direction.
diffx = vecnorm(diff(pc.Location,1,1), 2, 3);
diffy = vecnorm(diff(pc.Location,1,2), 2, 3);
sdx = size(diffx);
sdy = size(diffy);

% Determine the index addends corresponding to each configuration in the
% Cartesian distance matrices.
[dnxx,dnxy] = cellfun(@(ci) find(diff(ci,1,1)==0 & ci(1:2,:)), c, ...
    'UniformOutput', false);
[dnyx,dnyy] = cellfun(@(ci) find(diff(ci,1,2)==0 & ci(:,1:2)), c, ...
    'UniformOutput', false);

% Loop over all locations and configurations and invalidate all plane seeds
% that contain points outside the map or that contain neighboring points
% that are more than the specified maximum Cartesian distance apart.
parfor ip = 1 : np
    % Determine the subscript index of the point.
    [ox,oy] = ind2sub(spc, ip);
    
    % Loop over all configurations.
    for ic = 1 : nc
        % Determine the points that make up the plane seed.
        ocx = ox + dcx(ic,:); %#ok<PFBNS>
        ocy = oy + dcy(ic,:); %#ok<PFBNS>
        
        % Check if all points are located inside the map and if the
        % Cartesian distances between the point pairs fall below the
        % maximum admissible distance.
        if all([ocx,ocy]>=1 & [ocx,ocy]<=repelem(spc,1,4))
            % Determine the vertical and horizontal Cartesian distances.
            dx = diffx(sub2ind(sdx, ox+dnxx{ic}-2, oy+dnxy{ic}-2)); ...
                %#ok<PFBNS>
            dy = diffy(sub2ind(sdy, ox+dnyx{ic}-2, oy+dnyy{ic}-2)); ...
                %#ok<PFBNS>
            
            % Invalidate all point pairs whose Cartesian distances exceed
            % the specified maximum.
            if any([dx; dy] > lmax)
                % Invalidate the configuration.
                crte(ip,ic) = Inf;
            else
                % Store the linear point indices.
                crto{ip,ic} = uint32(sub2ind(spc, ocx, ocy));
            end
        else
            % Invalidate the configuration.
            crte(ip,ic) = Inf;
        end
    end
end

% Set the optimizer parameters.
cpuopt = optimoptions('lsqcurvefit', ...
    'Algorithm', 'levenberg-marquardt', ...
    'SpecifyObjectiveGradient', true, 'Display', 'none');
gpuopt.tolerance = 1e-6;
gpuopt.iterations = 20;
gpuopt.estimator = 0; % Least-squares estimator.
gpuopt.model = 8; % Raytracing model.

% Compute the error corresponding to each plane seed.
if strcmpi(device, 'mat') || strcmpi(device, 'cpu') % Optimization on CPU.
    % For each location and each configuration, compute the error of
    % creating the corresponding plane.
    parfor ip = 1 : np
        for ic = 1 : nc
            if isnan(crte(ip,ic))
                % Determine the initial plane parameters.
                x0 = plnx(ip,:);

                % Determine the input data required to compute the
                % intersections between the rays and the plane, i.e.
                % the ray direction vectors.
                xdata = v([s(ip,:),crto{ip,ic}],:); %#ok<PFBNS>

                % Retrieve the corresponding measured ray radii.
                ydata = r(crto{ip,ic}); %#ok<PFBNS>

                % Compute the error that the creation of the plane
                % would incur.
                if strcmpi(device, 'mat') % Native MATLAB optimization
                    [crtx{ip,ic},crte(ip,ic)] = lsqcurvefit(@rxp, ...
                        x0, xdata, ydata, [], [], cpuopt);
                else % Cpufit optimization
                    [parameters,~,crte(ip,ic)] = ...
                        cpufit(ydata, [], gpuopt.model, x0', ...
                        gpuopt.tolerance, gpuopt.iterations, [], ...
                        gpuopt.estimator, xdata', 8*numel(xdata));
                    crtx{ip,ic} = parameters';
                end                    
            end
        end
    end
else % Optimization on GPU.
    % Determine the initial plane parameters.
    [ip,ic] = find(isnan(crte));
    x0 = plnx(ip,:)';

    % Determine the input data required to compute the intersections
    % between the rays and the plane, i.e. the ray direction vectors.
    icrt = sub2ind(size(crto),ip,ic);
    o = vertcat(crto{icrt});
    xdata = reshape(v([s(ip,:),o]',:)', (3+4)*3, []);

    % Retrieve the corresponding measured ray radii.
    ydata = r(o');

    % Compute the error that the creation of the plane would incur.
    [parameters,~,crte(icrt)] ...
        = gpufit(ydata, [], gpuopt.model, x0, gpuopt.tolerance, ...
            gpuopt.iterations, [], gpuopt.estimator, xdata);
    crtx(icrt) = deal(cellfun(@transpose, num2cell(parameters,1), ...
        'UniformOutput', false));
    clear ic ip x0 o xdata ydata parameters
end

% For each point in the plane array, compute the indices of the plane seeds
% it is part of.
for icrt = 1 : numel(crto)
    for oi = crto{icrt}
        plncrt{oi}(end+1) = icrt;
    end
end

%% Create, extend, and merge planes.
% Initialize the map.
m = NaN(spc);

% Create a colormap of different colors.
colormap = hsv(np);
rng(0)
colormap = colormap(randperm(np),:);

% Plot the point cloud.
if any(strcmpi(display, {'iter','final'}))
    % Plot the point cloud.
    cdata = ones(size(l)) / 2;
    sctr = scatter3(l(:,1), l(:,2), l(:,3), sqrt(np)/10, cdata, '.');

    % Configure the figure.
    fig = sctr.Parent.Parent;
    fig.Name = mfilename;
    if strcmpi(display, 'iter')
        fig.Visible = 'on';
    else
        fig.Visible = 'off';
    end
    labelaxes
    axis equal
    campos([0,0,0])
end

% Initialize the set of potential plane extensions.
exte = NaN(0,1);
extelb = NaN(0,1);
exto = NaN(0,1);
extp = NaN(0,1);
extx = NaN(0,3);

% Initialize the set of potential mergers of planes.
mrge = NaN(0,1);
mrgelb = NaN(0,1);
mrgp = NaN(0,2);
mrgx = NaN(0,3);

% Determine the action that incurs minimum incremental error.
[demin,icrtemin] = min(crte(:));
action = 1;

% Initialize the cell array where to store the plane extraction results.
res = repmat(struct('plane', [], 'steps', [], 'emax', []), 0, 1);

% Remove all infeasible stopping criteria.
nstop = nstop(nstop<np);
estop = estop(estop>=demin);

% Initialize the number of planes.
npln = np;

% Create, extend, merge planes and compute the error corresponding to the
% next step until the least restrictive stopping criterion is met.
it = 1;
while isfinite(demin) && ((~isempty(estop) && demin <= estop(end)) ...
        || (~isempty(nstop) && npln > nstop(end)))    
    % Perform the action that incurs the smallest incremental error.
    if action == 1 % Create plane.
        % Update the number of planes.
        npln = npln - 3;
        
        % Determine which plane to create, the error incurred by creating
        % the plane, the points the plane consists of, and the plane
        % parameters.
        [p,~] = ind2sub(size(crte), icrtemin);
        e = crte(icrtemin);
        o = crto{icrtemin};
        x = crtx{icrtemin};
    elseif action == 2 % Extend plane.
        % Update the number of planes.
        npln = npln - 1;
    
        % Determine which plane to extend, the error corresponding to the
        % extended plane, the point added to the plane, and the parameters
        % of the extended plane.
        p = extp(iextdemin);
        e = exte(iextdemin);
        o = exto(iextdemin);
        x = extx(iextdemin,:);
    end

    % If a plane was created or extended, update the corresponding data.
    % Otherwise, merge planes.
    if action <= 2 % Update data.
        % Update the map.
        m(o) = p;
        
        % Update the created or extended plane.
        plne(p) = e;
        plno{p} = [plno{p},o];
        plnx(p,:) = x;
         
        % Invalidate all superseded plane seeds.
        crte([plncrt{o}]) = Inf;
        
        % Remove all superseded plane extensions.
        iext = ~ismember(exto,o);
        exte = exte(iext,:);
        extelb = extelb(iext,:);
        exto = exto(iext,:);
        extp = extp(iext,:);
        extx = extx(iext,:);
        
        % Update all potential extensions of the plane.
        iext = extp==p;
        exte(iext,:) = NaN;
        extelb(iext,:) = max(max(exte(iext),extelb(iext)), plne(p));
        extx(iext,:) = repmat(plnx(p,:), sum(iext), 1);
        
        % Determine the indices of all points in the 4-neighborhood of the
        % added points. Remember how they extend the plane, e.g. in
        % positive x-direction, in negative x-direction, in positive
        % y-direction, etc.
        [ox,oy] = ind2sub(spc, o');
        ox = ox + [-1,1,0,0];
        oy = oy + [0,0,-1,1];
        ext = repmat(1:4, size(o'));
        
        % Remove all neighboring points that are located outside the map.
        io = ox>=1 & ox<=spc(1) & oy>=1 & oy<=spc(2);
        ox = ox(io);
        oy = oy(io);
        ext = ext(io);
        
        % Remove all neighboring points that do not correspond to
        % reflections.
        io = isfinite(r(sub2ind(spc, ox, oy)));
        ox = ox(io);
        oy = oy(io);
        ext = ext(io);
                
        % Find out which extensions exceed the maximum Cartesian distance.
        oxxneg = ox(ext==1);
        oyxneg = oy(ext==1);
        ioxneg = diffx(sub2ind(sdx,oxxneg,oyxneg)) > lmax;
        oxneg = sub2ind(spc, oxxneg(ioxneg), oyxneg(ioxneg));
        oxxpos = ox(ext==2);
        oyxpos = oy(ext==2);
        ioxpos = diffx(sub2ind(sdx,oxxpos-1,oyxpos)) > lmax;
        oxpos = sub2ind(spc, oxxpos(ioxpos), oyxpos(ioxpos));
        oxyneg = ox(ext==3);
        oyyneg = oy(ext==3);
        ioyneg = diffy(sub2ind(sdy,oxyneg,oyyneg)) > lmax;
        oyneg = sub2ind(spc, oxyneg(ioyneg), oyyneg(ioyneg));
        oxypos = ox(ext==4);
        oyypos = oy(ext==4);
        ioypos = diffy(sub2ind(sdy,oxypos,oyypos-1)) > lmax;
        oypos = sub2ind(spc, oxypos(ioypos), oyypos(ioypos));
        
        % Remove all neighboring points that exceed the maximum Cartesian
        % distance.
        on = unique(sub2ind(spc, ox, oy));
        on = on(~ismember(on, ...
            unique([oxneg(:); oxpos(:); oyneg(:); oypos(:)])));
        
        % Determine which neighboring points represent new extensions.
        onext = on(~ismember(on,exto(iext)) & isnan(m(on)));
        
        % Add the new extensions.
        iext = numel(exte) + (1:numel(onext));
        exte(iext,:) = NaN;
        extelb(iext,:) = e;
        exto(iext,:) = onext;
        extp(iext,:) = p;
        extx(iext,:) = repmat(x, numel(iext), 1);

        % Determine the neighboring planes.
        pn = unique(m(on(isfinite(m(on)) & m(on)~=p)));
        
        % Estimate the errors incurred by merging the plane with each of
        % the neighboring planes.
        for ipn = 1 : numel(pn)
            % Check if a merger for the plane pair already exists.
            pm = sort([p,pn(ipn)]);
            imrg = find(all(mrgp==pm, 2));
            
            % If a merger does not yet exist, add it. Otherwise, update it.
            if isempty(imrg)
                imrg = numel(mrge) + 1;
                mrgp(imrg,:) = pm; %#ok<*AGROW>
                mrgelb(imrg,:) = sum(plne(pn));
                mrgx(imrg,:) = plnx(p,:);
            else
                mrgelb(imrg,:) ...
                    = max([mrge(imrg),mrgelb(imrg),sum(plne(pm))]);
            end
            mrge(imrg,:) = NaN;
        end
    else % Merge planes.
        % Update the number of planes.
        npln = npln - 1;
    
        % Update the map.
        pm = mrgp(imrgdemin,:);
        m(m==pm(2)) = pm(1);
        
        % Find all planes that are neighbors of either of the merged
        % planes.
        pn = unique(mrgp(any(ismember(mrgp,pm),2),:));
        pn = pn(~ismember(pn,pm));
        
        % Loop over all neighboring planes and update the corresponding
        % mergers.
        for ipn = 1 : numel(pn)
            % Determine the indices of the mergers that connect the first
            % plane and the second plane to the neighboring plane.
            imrg1n = all(mrgp==sort([pm(1),pn(ipn)]), 2);
            imrg2n = all(mrgp==sort([pm(2),pn(ipn)]), 2);
            
            % Compute the lower bound of the merging error.
            elb = max([mrge(imrgdemin) + plne(pn(ipn)); ...
                [mrge(imrg1n); mrgelb(imrg1n)] + plne(pm(2)); ...
                [mrge(imrg2n); mrgelb(imrg2n)] + plne(pm(1))]);
            
            % Update or add the merger.
            if all(~imrg1n)
                imrg1n = numel(mrge) + 1;
                mrgp(imrg1n,:) = sort([pm(1),pn(ipn)]);
                mrgx(imrg1n,:) = mrgx(imrgdemin,:);
            end
            mrge(imrg1n,:) = NaN;
            mrgelb(imrg1n,:) = elb;
        end
        
        % Determine the potential extensions of the merged plane.
        iext1 = find(extp==pm(1));
        iext2 = find(extp==pm(2));
        oext = unique(exto([iext1; iext2]));
        [~,iexto1] = ismember(oext, exto(iext1));
        [~,iexto2] = ismember(oext, exto(iext2));
        
        % Compute the lower bounds of the errors corresponding to extending
        % the merged plane.
        elb1 = zeros(size(oext));
        elb2 = zeros(size(oext));
        elb1(iexto1>0) = max(extelb(iext1(iexto1(iexto1>0))), ...
            exte(iext1(iexto1(iexto1>0)))) + plne(pm(2));
        elb2(iexto2>0) = max(extelb(iext2(iexto2(iexto2>0))), ...
            exte(iext2(iexto2(iexto2>0)))) + plne(pm(1));
        
        % Update the errors and plane parameters corresponding to extending
        % the merged plane.
        iext = [iext1; iext2];
        exte(iext,:) = [];
        extelb(iext,:) = [];
        exto(iext,:) = [];
        extp(iext,:) = [];
        extx(iext,:) = [];
        iext = numel(exte) + (1:numel(oext));
        exte(iext,:) = NaN;
        extelb(iext,:) = max(elb1,elb2);
        exto(iext,:) = oext;
        extp(iext,:) = pm(1);
        extx(iext,:) = repmat(mrgx(imrgdemin,:), numel(iext), 1);
        
        % Update the plane array.
        plne(pm(1),:) = mrge(imrgdemin);
        plno{pm(1),:} = [plno{pm}];
        plnx(pm(1),:) = mrgx(imrgdemin,:);
        plne(pm(2),:) = NaN;
        
        % Remove all mergers for the second plane.
        imrg = ~any(mrgp==pm(2), 2);
        mrge = mrge(imrg,:);
        mrgelb = mrgelb(imrg,:);
        mrgp = mrgp(imrg,:);
        mrgx = mrgx(imrg,:);
    end
    
    % Find out which plane in which configuration yields minimum error when
    % being created.
    [crtemin,icrtemin] = min(crte(:));
    crtemin(~isfinite(crtemin)) = Inf;
    
    % Compute the error increments corresponding to extending the planes
    % and find their minimum.
    extde = exte - plne(extp);
    extdelb = extelb - plne(extp);
    [extdemin,iextdemin] = min(extde);
    [extdelbmin,iextdelbmin] = min(extdelb);
    extdemin(isempty(extdemin)) = Inf;
    extdelbmin(isempty(extdelbmin)) = Inf;
    
    % Compute the error increments corresponding to merging planes and find
    % their minimum.
    em = plne(mrgp(:,1)) + plne(mrgp(:,2));
    mrgde = mrge - em;
    mrgdelb = mrgelb - em;
    [mrgdemin,imrgdemin] = min(mrgde);
    [mrgdelbmin,imrgdelbmin] = min(mrgdelb);
    mrgdemin(isempty(mrgdemin)) = Inf;
    mrgdelbmin(isempty(mrgdelbmin)) = Inf;
    
    % Recompute all extension errors and merging errors that could lead to
    % smaller error increments than the minimum error increment determined
    % so far.
    while min([extdelbmin,mrgdelbmin]) < min([crtemin,extdemin,mrgdemin])
        % Recompute the smallest lower bound incremental error.
        if extdelbmin < mrgdelbmin % Recompute extension error.
            % Determine the indices of the points represented by the
            % extended plane.
            pext = extp(iextdelbmin);
            o = [plno{pext},exto(iextdelbmin)];
            
            % Determine the initial plane parameters.
            x0 = extx(iextdelbmin,:);

            % Determine the input data required to compute the
            % intersections between the rays and the plane, i.e. the ray
            % direction vectors.
            xdata = v([s(pext,:),o],:);

            % Retrieve the corresponding measured ray radii.
            ydata = r(o);

            % Compute the exact error that the extension of the plane
            % incurs.
            if strcmpi(device, 'mat')
                [extx(iextdelbmin,:),exte(iextdelbmin)] ...
                    = lsqcurvefit(@rxp, x0, xdata, ydata, [], [], cpuopt);
            else
                [parameters,~,exte(iextdelbmin)] = ...
                    cpufit(ydata, [], gpuopt.model, x0', ...
                    gpuopt.tolerance, gpuopt.iterations, [], ...
                    gpuopt.estimator, xdata', 8*numel(xdata));
                extx(iextdelbmin,:) = parameters';
            end

            % Update the lower bound of the extension error.
            extelb(iextdelbmin) = exte(iextdelbmin);
            
            % Update the error increments.
            extde(iextdelbmin) = exte(iextdelbmin) - plne(pext);
            extdelb(iextdelbmin) = extelb(iextdelbmin) - plne(pext);
            [extdemin,iextdemin] = min(extde);
            [extdelbmin,iextdelbmin] = min(extdelb);
        else % Recompute merging error.
            % Determine the indices of the planes.
            pm = mrgp(imrgdelbmin,:);
            
            % Determine the indices of the points represented by either
            % plane.
            o = [plno{pm}];
            
            % Determine the initial plane parameters.
            x0 = mrgx(imrgdelbmin,:);
            
            % Determine the input data required to compute the
            % intersections between the rays and the plane, i.e. the ray
            % direction vectors.
            xdata = v([s(pm(1),:),o],:);
            
            % Retrieve the corresponding measured ray radii.
            ydata = r(o);
            
            % Compute the exact merging error.
            if strcmpi(device, 'mat')
                [mrgx(imrgdelbmin,:),mrge(imrgdelbmin)] ...
                    = lsqcurvefit(@rxp, x0, xdata, ydata, [], [], cpuopt);
            else
                [parameters, ~, mrge(imrgdelbmin)] = ...
                    cpufit(ydata, [], gpuopt.model, x0', ...
                    gpuopt.tolerance, gpuopt.iterations, [], ...
                    gpuopt.estimator, xdata', 8*numel(xdata));
                mrgx(imrgdelbmin,:) = parameters';
            end
            % Update the lower bound of the merging error.
            mrgelb(imrgdelbmin) = mrge(imrgdelbmin);
            
            % Update the error increments.
            em = sum(plne(pm));
            mrgdelb(imrgdelbmin) = mrgelb(imrgdelbmin) - em;
            mrgde(imrgdelbmin) = mrge(imrgdelbmin) - em;
            [mrgdelbmin,imrgdelbmin] = min(mrgdelb);
            [mrgdemin,imrgdemin] = min(mrgde);
        end
    end
    
    % Find the action that incurs the smallest incremental error.
    [demin,action] = min([crtemin,extdemin,mrgdemin]);
    
    % If the error increment corresponding to the next step is higher than
    % the current error threshold or if the plane count fell below a count
    % specified in the stopping criteria, store the current plane
    % extraction result.
    if ~isfinite(demin) || (~isempty(estop) && demin>estop(1)) ...
            || (~isempty(nstop) && npln<=nstop(1))
        % Create cell arrays that contain the information of the regular
        % planes extracted so far.
        [~,irpln] = sort(plne);
        irpln = irpln(isfinite(plne(irpln)));
        rplnec = num2cell(plne(irpln));
        rplnxc = mat2cell(v(reshape(s(irpln,:)',[],1),:) ...
            .* reshape(plnx(irpln,:)',[],1), 3*ones(numel(irpln),1));
        rplnoc = plno(irpln);
        
        % Create cell arrays that contain the information of the atomic
        % planes extracted so far.
        iapln = find(isnan(m));
        aplnec = num2cell(zeros(size(iapln)));
        aplnxc = mat2cell(v(reshape(s(iapln,:)',[],1),:) ...
            .* reshape(plnx(iapln,:)',[],1), 3*ones(numel(iapln),1)); 
        aplnoc = num2cell(iapln);
        
        % Create a struct that stores the planes extracted so far and add
        % it to the results.
        pln = repmat(struct('index',[],'error',[],'param',[]), ...
            numel(irpln)+numel(iapln), 1);
        [pln.index] = deal(rplnoc{:}, aplnoc{:});
        [pln.error] = deal(rplnec{:}, aplnec{:});
        [pln.param] = deal(rplnxc{:}, aplnxc{:});
        emax = [];
        if ~isempty(estop)
            emax = estop(1);
        end
        res(end+1) = struct('plane', pln, 'emax', emax, 'steps', it);
        
        % Remove outdated stopping criteria.
        estop(estop<demin) = [];
        nstop(nstop>=npln) = [];
    end
    
    % Update the plot.
    if strcmpi(display, 'iter')
        sctr.CData(isfinite(m(:)),:) = colormap(m(isfinite(m(:))),:);
        drawnow limitrate
    end
    
    % Check data consistency.
    %check(m, plno, plne, plnx, crte, crto, crtx, ...
    %    exte, extelb, exto, extp, extx, mrge, mrgelb, mrgp, mrgx)
    
    % Increment the action iterator.
    it = it + 1;
end

%% Visualize result.
% Visualize the points and planes.
if any(strcmpi(display, {'iter','final'}))
    % Show the point cloud.
    sctr.CData(isfinite(m(:)),:) = colormap(m(isfinite(m(:))),:);
    fig.Visible = 'on';
    
    % Plot the planes.
    h = ishold;
    hold(sctr.Parent, 'on')
    for ip = find(isfinite(plne))'
        % Determine the subscript indices of the points of the plane.
        [ox,oy] = ind2sub(size(m), plno{ip});
        
        % Determine the indices of the points that form the boundary of the
        % plane.
        k = boundary(ox(:), oy(:));
        
        % Determine the 3-D coordinates of the boundary.
        lk = v(plno{ip}(k),:) ...
            .* rxp(plnx(ip,:), v([s(ip,:),plno{ip}(k)],:));
        
        % Plot the boundary.
        fill3(lk(:,1), lk(:,2), lk(:,3), colormap(ip,:), 'FaceAlpha', 0.8)
    end
    if ~h
        hold(sctr.Parent, 'off')
    end
end

end

function ind = sub2ind(sz, i, j)
% SUB2IND Linear index from two subscripts.
%   SUB2IND is an accelerated implementation of MATLAB's SUB2IND function
%   for two subscript indices.

ind = (j-1)*sz(1) + i;

end

function check(m, plno, plne, plnx, crte, crto, crtx, ...
    exte, extelb, exto, extp, extx, mrge, mrgelb, mrgp, mrgx) %#ok<DEFNU>
% CHECK Check consistency of given matrices.
%   CHECK makes sure the given matrices are consistent. It is meant for
%   debugging and testing the function PCEXTRPLN.

% Check the map.
assert(all(mod(m(isfinite(m)),1)==0))
assert(all(isnan(m(~isfinite(m)))))
assert(all(m(isfinite(m))>=1))
assert(all(m(isfinite(m))<=numel(m)))

% Check the planes.
assert(isequal(numel(plno),numel(plne),size(plnx,1)))
assert(size(plnx,2)==3)

% Check consistency of map and planes.
assert(numel(unique(m(isfinite(m))))==sum(isfinite(plne)))
assert(numel(m(isfinite(m)))==sum(cellfun(@numel,plno(isfinite(plne)))))
for ip = 1 : numel(m)
    if isfinite(plne(ip))
        assert(all(m(plno{ip})==ip))
    end
end

% Check the plane seeds.
assert(isequal(size(crte),size(crto),size(crtx)))
assert(all(crte(:)>=0))

% Check the plane extensions.
assert(isequal(size(exte),size(extelb),size(exto),size(extp)))
assert(size(extx,1)==numel(exte))
assert(size(extx,2)==3)

% Check consistency of map and plane extensions.
assert(all(ismember(extp,m)))
assert(all(exto>=1 & exto<=numel(m)))

% Check the plane mergers.
assert(isequal(size(mrge),size(mrgelb)))
assert(isequal(size(mrgp,1),size(mrgx,1),numel(mrge)))
assert(size(mrgp,2)==2)
assert(size(mrgx,2)==3)

% Check consistency of map and mergers.
assert(all(ismember(mrgp(:),m)))

end
