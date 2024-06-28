function wf = spm_pac_transform_volume(f, y, itrp, bnd)
% FORMAT wf = transform_volume(f, y, [itrp], [bnd])
%   f    - Input volume [Nx Ny Nz ...]
%   y    - Transformation [Mx My Mz 3]
%   itrp - Interpolation order {default: 1}
%   bnd  - Boundary conditions {default: 0 = no wrapping around}
%   wf   - Warped volume [Mx My Mz ...]
%
% Note that the transformation must map two voxel spaces, so it should have
% already been composed with voxel-to-world affine matrices.
%% Yael Balbastre 2024

if nargin < 4
    bnd = 0;
    if nargin < 3
        itrp = 1;
    end
end

dimin  = size(f);
f      = reshape(f, dimin(1), dimin(2), dimin(3), []);
dimout = size(y);

wf = zeros([dimout(1:3) size(f,4)], 'like', f);
for k=1:size(f,4)
    c = spm_diffeo('bsplinc', single(f(:,:,:,k)), [itrp itrp itrp bnd bnd bnd]);
    %wf(:,:,:,k) = spm_diffeo('bsplins', c, single(numeric(y)), [itrp itrp itrp bnd bnd bnd]);	% JD: Odd error with numeric
	wf(:,:,:,k) = spm_diffeo('bsplins', c, single(full(y)), [itrp itrp itrp bnd bnd bnd]);
end

wf = reshape(wf, [dimout(1:3) dimin(4:end)]);