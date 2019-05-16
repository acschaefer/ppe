function pln = extrplnmsac(pc,x)
% EXTRPLNMSAC Extract planes from point cloud using MSAC algorithm.
%   IDX = EXTRPLNMSAC(PC,X) extracts planes from the point cloud PC using
%   the M-estimator sample consensus algorithm.
%
%   X is a 2-element vector. The first element specifies the maximum
%   distance from an inlier point to the plane, the second controls the
%   percentage of remaining unassigned points after plane extration.
%
%   PLN is a struct array. PLN(m).param tells the parameters of the m-th
%   plane, PLN(m).index tells the indices that represent the plane.
%
%   Example:
%      load('object3d.mat')
%      pln = extrplnmsac(ptCloud, [0.05, 0.10])
%
%   See also PCFITPLANE, PCEXTRPLN.

% Copyright 2018 Alexander Schaefer

% Rearrage point cloud.
pc = pointCloud(reshape(pc.Location,[],3,1));

% Extract planes from point cloud until a specified percentage of points is
% assigned to planes.
pln = [];
while sum(all(isfinite(pc.Location),2)) > max(x(2)*pc.Count,2)
    % Fit a plane to the points in the point cloud.
    [model,in] = pcfitplane(pc, x(1));
    
    % Invalidate the points that belong to the plane.
    p = pc.Location;
    p(in,:) = NaN;
    pc = pointCloud(p);
    
    % Add the indices of the points and plane parameters to the result.
    pln(end+1).index = in; %#ok<AGROW>
    o = [-model.Parameters(4)/model.Parameters(1),0,0];
    r = [-model.Parameters(2)/model.Parameters(1),1,0];
    s = [-model.Parameters(3)/model.Parameters(1),0,1];
    pln(end).param = [o; o+r; o+s];
end

end
