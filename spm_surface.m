function spm_surface(P,Q,U,V,B)
% renders SPM{Z} on a surface
% FORMAT spm_surface(P,Q,U,V,B)
% P  - filename of 8 bit SPM{Z} image
% Q  - filename of 8 bit image to surface render
% U  - threshold for SPM{Z}
% V  - threshold for reference image {0 - 255}
% B  - transformation matrix
%____________________________________________________________________________
%
% spm_surface renders surprathreshold regions from a SPM{Z} in image [ANALYZE]
% format (created by spm_spm.m) beneath the cortical surface of a specified
% reference image.
%
% The reference image is usually a structural MRI image and should be
% 8 bit {unsigned char} data.  Both images require properly specifed
% ORIGINS in the header (this is guaranteed if the data have been spatially
% normalized)
%
% If the original data conform to a radiological convention the results
% of this rendering will appear mirror reversed.
%
%__________________________________________________________________________
% %W% %E%

%----------------------------------------------------------------------------
global CWD

% memory map reference image and create transformation matrix {A}
%----------------------------------------------------------------------------
[d d SCALEZ d d ORIGINZ] = spm_hread(P);			% get origin
[d d SCALES d d ORIGINS] = spm_hread(Q);			% get origin

B0     = spm_matrix([0 0 0 0 pi/2 pi]);			% default view
U      = U/SCALEZ;					% Z threshold
M      = 128;						% display -
N      = 128;						% dimensions {mm}
Vz     = spm_map(P);					% memory map
Vs     = spm_map(Q);					% memory map
Az     = spm_matrix(-ORIGINZ);				% center on origin
As     = spm_matrix(-ORIGINS);				% center on origin
Az     = spm_matrix([0 0 0 0 0 0 Vz(4:6)'/2])*Az;	% anisotropy
As     = spm_matrix([0 0 0 0 0 0 Vs(4:6)'/2])*As;	% anisotropy
Az     = B*B0*Az;					% implement transform
As     = B*B0*As;					% implement transform
Az     = spm_matrix([M/2 N/2 M/2])*Az;			% re-center to corner
As     = spm_matrix([M/2 N/2 M/2])*As;			% re-center to corner


% get surface rendering and indices
%----------------------------------------------------------------------------
RENZ = spm_render_vol(Vz,Az,[M N],[U 2]);
RENS = spm_render_vol(Vs,As,[M N],[V 2]);

% normalize and combine %----------------------------------------------------------------------------
RENZ   = 65*(RENZ > 0) + 48*RENZ/max(RENZ(:));
RENS   = 64*(RENS/max(RENS(:)) + (~RENS));
REN    = spm_merge(RENZ,RENS);

% display %----------------------------------------------------------------------------
figure(3); delete(gca)
load Split; colormap(split)
image(REN); axis xy; axis image; axis off
title('Surface rendered SPM{Z}')
text(24,-16,'results directory:')
text(96,-16,CWD,'Fontsize',16,'Fontweight','Bold')

% unmap %----------------------------------------------------------------------------
spm_unmap_vol(Vz)
spm_unmap_vol(Vs)
