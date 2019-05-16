function labelaxes(unit)
% LABELAXES Label x, y, and z axis.
%   LABELAXES() labels the x, y, and z axis with 'x', 'y', and 'z'. 
%
%   LABELAXES(UNIT) adds the UNIT string to each label. For example, if
%   UNIT is 'm', the x-label is 'x [m]'.

% Copyright 2016-2017 Alexander Schaefer

%% Validate input.
% Define the unit, if not given.
if nargin < 1
    unit = '';
end

% Check the input argument.
if ~ischar(unit)
    error('UNIT must be a string.')
end

%% Write labels.
if isempty(unit)
    xlabel('x')
    ylabel('y')
    zlabel('z')
else
    xlabel(['x [', unit, ']'])
    ylabel(['y [', unit, ']'])
    zlabel(['z [', unit, ']'])
end

end
