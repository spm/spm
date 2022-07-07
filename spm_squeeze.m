function B = spm_squeeze(A, dim)
% Version of squeeze with the possibility to select the dimensions to remove
% FORMAT  B = spm_squeeze(A, dim)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2010-2022 Wellcome Centre for Human Neuroimaging


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
