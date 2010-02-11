function B = spm_squeeze(A, dim)
% version of squeeze with the possibility to select the dimensions to remove
% FORMAT  B = spm_squeeze(A, dim)
%
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_squeeze.m 3722 2010-02-11 16:23:28Z vladimir $

if nargin == 1
    B = squeeze(A);
else
    siz = size(A);
    dim = intersect(dim, find(siz == 1));
    if ~isempty(dim)
        siz(dim) = [];
        B = reshape(A, siz);
    else
        B = A;
    end
end
