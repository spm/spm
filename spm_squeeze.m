function B = spm_squeeze(A, dim)
% version of squeeze with the possibility to select the dimensions to remove
% FORMAT  B = spm_squeeze(A, dim)
%
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_squeeze.m 5794 2013-12-09 12:41:52Z vladimir $

if nargin == 1
    B = squeeze(A);
else
    siz = size(A);
    dim = intersect(dim, find(siz == 1));
    if ~isempty(dim)
        siz(dim) = [];
        if size(siz) == 1
            siz = [siz 1];
        end
        B = reshape(A, siz);
    else
        B = A;
    end
end
