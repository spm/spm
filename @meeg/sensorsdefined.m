function res = sensorsdefined(this, type)
% Returns 1 if there are sensors defined in the obgect, 0 otherwise
% FORMAT res = sensorsdefined(this, type)
%   type (optional) - 'EEG' or 'MEG'
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: sensorsdefined.m 1305 2008-04-03 17:39:04Z vladimir $

if isempty(this.sensors) || isempty(this.sensors.type)
    res = 0;
    return;
else
    if nargin == 1
        res = 1;
        return;
    else
        res = ~isempty(strmatch(type, this.sensors.type, 'exact'));
    end
end
    
