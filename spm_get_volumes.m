function gl = spm_get_volumes(P)
% Compute total volumes from grey matter images
% FORMAT gl = spm_get_volumes(P)
% gl - a vector of volumes (in litres)
% P  - a matrix of image filenames
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_get_volumes.m 1131 2008-02-06 11:17:09Z spm $


if nargin<1,
    P = spm_select(Inf,'image','Select images');
end;
V = spm_vol(P);
n = numel(V);
gl = zeros(n,1);
spm_progress_bar('Init',n,'Computing volumes');
for j=1:n,
    tot = 0;
    for i=1:V(j).dim(3),
        img = spm_slice_vol(V(j),spm_matrix([0 0 i]),V(j).dim(1:2),0);
        img = img(isfinite(img));
        tot = tot + sum(img(:));
    end;
    gl(j) = tot*abs(det(V(j).mat)/100^3);
    spm_progress_bar('Set',j);
end;
spm_progress_bar('Clear');
