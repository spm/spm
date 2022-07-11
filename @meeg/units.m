function res = units(this, varargin)
% Method for setting/getting all units, over channels
% FORMAT res = units(this, ind)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


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
