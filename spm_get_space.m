function M = spm_get_space(P,M)
% Get/set the voxel-to-world mapping of an image
% FORMAT M = spm_get_space(P)
%            spm_get_space(P,M)
% M - voxel-to-world mapping
% P - image filename
%_______________________________________________________________________
% John Ashburner $Id$

[pth,nam,ext] = fileparts(P);
t = find(ext==',');
n = [1 1];
if ~isempty(t),
    n   = str2num(ext((t(1)+1):end));
    ext = ext(1:(t-1));
    P   = fullfile(pth,[nam ext]);
end;

N = nifti(P);
if nargin==2,
    if n(1)==1,
        N.mat = M;
        if ~isempty(N.extras) && isstruct(N.extras) && isfield(N.extras,'mat') &&...
            size(N.extras.mat,3)>=1 && ~sum(N.extras.mat(:,:,1)),
            N.extras.mat(:,:,1) = M;
        end;
    else
        if sum((N.mat(:)-M(:)).^2) > 1e-4,
            N.extras.mat(:,:,n(1)) = M;
        end;
    end;
    create(N);
else
    if ~isempty(N.extras) && isstruct(N.extras) && isfield(N.extras,'mat') &&...
        size(N.extras.mat,3)>=n(1) && sum(sum(N.extras.mat(:,:,n(1)))),
        M = N.extras.mat(:,:,n(1));
    else
        M = N.mat;
    end;
end;
