function res = units(this, varargin)
% Method for setting/getting all units, over channels
% FORMAT res = units(this, ind)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: units.m 4432 2011-08-15 12:43:44Z christophe $

if this.montage.Mind == 0
    res = getset(this, 'channels', 'units', varargin{:});
else
    res = getset(this.montage.M(this.montage.Mind), 'channels', 'units', varargin{:});
end