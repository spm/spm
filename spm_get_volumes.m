function gl = spm_get_volumes(P)
% Compute total volumes from tissue segmentations
% FORMAT gl = spm_get_volumes(P)
% P  - a matrix of image filenames
% gl - a vector of volumes (in litres)
%__________________________________________________________________________
% Copyright (C) 2006-2011 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_get_volumes.m 4352 2011-06-13 17:27:46Z ged $


if ~nargin
    [P,sts] = spm_select(Inf,'image','Select images');
    if ~sts, gl = []; return; end
end
if iscellstr(P), P = char(P); end
V  = spm_vol(P);
n  = numel(V);
gl = zeros(n,1);
spm_progress_bar('Init',n,'Computing volumes');
for j=1:n
    tot = 0;
    for i=1:V(j).dim(3)
        img = spm_slice_vol(V(j),spm_matrix([0 0 i]),V(j).dim(1:2),0);
        img = img(isfinite(img));
        tot = tot + sum(img(:));
    end
    gl(j) = tot*abs(det(V(j).mat)/100^3);
    spm_progress_bar('Set',j);
end
spm_progress_bar('Clear');
