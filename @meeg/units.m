function res = units(this, varargin)
% Method for setting/getting all units, over channels
% FORMAT res = units(this, ind)
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: units.m 5025 2012-10-31 14:44:13Z vladimir $

if this.montage.Mind == 0
    res = getset(this, 'channels', 'units', varargin{:});
else
    res = getset(this.montage.M(this.montage.Mind), 'channels', 'units', varargin{:});
end