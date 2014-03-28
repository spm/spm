function res = units(this, varargin)
% Method for setting/getting all units, over channels
% FORMAT res = units(this, ind)
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: units.m 5933 2014-03-28 13:22:28Z vladimir $

if this.montage.Mind == 0
    res = getset(this, 'channels', 'units', varargin{:});
else
    if nargin == 3
        this.montage.M(this.montage.Mind) = getset(this.montage.M(this.montage.Mind), 'channels', 'units', varargin{:});
        res = this;
    else
        res = getset(this.montage.M(this.montage.Mind), 'channels', 'units', varargin{:});
    end
end