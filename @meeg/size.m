function res = size(this, varargin)
% returns the dimensions of the data matrix
% FORMAT res = size(this, dim))
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: size.m 4432 2011-08-15 12:43:44Z christophe $

if this.montage.Mind==0
    res = size(this.data.y, varargin{:});
else
    d_sz = size(this.data.y);
    d_sz(1) = size(this.montage.M(this.montage.Mind).tra,1);
    if isempty(varargin)
        res = d_sz;
    elseif varargin{:}<=length(d_sz)
        res = d_sz(varargin{:});
%     else
%         res = 1;
    end
end

if ntrials(this) == 1 && isempty(varargin)
    res = [res 1];
end
