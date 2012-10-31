function res = size(this, varargin)
% returns the dimensions of the data matrix
% FORMAT res = size(this, dim))
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: size.m 5025 2012-10-31 14:44:13Z vladimir $


if islinked(this)
    res = size(this.data, varargin{:});
    if this.montage.Mind~=0
        res(1) = size(this.montage.M(this.montage.Mind).tra,1);
    end
else
    if ~strncmpi(transformtype(this), 'TF', 2)
        res = zeros(1, 3);
    else
        res = zeros(1, 4);
    end
end


if ntrials(this) == 1 && isempty(varargin)
    res = [res 1];
end
