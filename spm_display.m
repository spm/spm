function spm_display(P,A)
% image display in three orthogonal sections
% FORMAT spm_display(P,[A])
% P    - filename
% A    - transformation matrix
%___________________________________________________________________________
%
% spm_display displays sagittal, coronal and transverse sections of the
% specified image {P}, through the center of the image space.  Sections
% are displayed in real space (accounting for anisotroptic voxel sizes)
% If a transformation matrix {A} is specified this transformation will
% be effected after centering the image and accounting for voxel anisotropy.
%
% The superimposed red lines are the planes that contain the ORIGIN as
% specified in the header.  To change this, or any other parameter, save
% the modified header (and display the image again to check any changes)
%
%__________________________________________________________________________
% %W% %E%


% if A is not specified set to the identity matrix
%----------------------------------------------------------------------------
if nargin < 2; A = spm_matrix([]); end

% display near the top of the screen H = [0 - 1]
%----------------------------------------------------------------------------
H      = 0.91;

% read header and memory map image
%----------------------------------------------------------------------------
[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(P);
V      = spm_map(P);


% apply A after centering the image and correcting for voxel anisotropy
%----------------------------------------------------------------------------
C      = spm_matrix(-(1+DIM)/2);				% image center
I      = spm_matrix([0 0 0 0 0 0 VOX]);			% isotropy
A      = inv(C)*inv(I)*A*I*C;				% transformation
ORIGIN = A*[ORIGIN(:); 1];


% read 3 orthogonal sections through the image center
%----------------------------------------------------------------------------
D      = [0 0 1 0;0 1 0 0;-1 0 0 DIM(1)/2;0 0 0 1]*A;
Ds     = spm_slice_vol(V,inv(D), DIM([3 2]),1);

D      = [0 0 1 0;1 0 0 0;0 -1 0 DIM(2)/2;0 0 0 1]*A;
Dc     = spm_slice_vol(V,inv(D), DIM([3 1]),1);

D      = [1 0 0 0;0 1 0 0;0 0 1 -(DIM(3) + 1)/2 ;0 0 0 1]*A;
Dt     = spm_slice_vol(V,inv(D), DIM([1 2]),1);

[DS I] = hist(Dt(:),32);

spm_unmap_vol(V);

% compute axes to correct for anisotropy of voxels and (normalized) window
%----------------------------------------------------------------------------
set(gcf,'Units','pixels')
WIN    = get(gcf,'Position');
WIN    = WIN(3)/WIN(4);
Y      = 0.36*DIM(2)*VOX(2)/max(DIM.*VOX);
X      = Y*DIM(1)*VOX(1)/(DIM(2)*VOX(2));
Z      = Y*DIM(3)*VOX(3)/(DIM(2)*VOX(2));


% display sections and draw lines though the center
%----------------------------------------------------------------------------
axes('Position',[0.1 (H - Z*WIN) Y Z*WIN])
imagesc(Ds); axis('xy')
title('sagittal')'; ylabel('z {voxels}')
line([0 DIM(2)],[1 1]*ORIGIN(3))
line([1 1]*ORIGIN(2),[0 DIM(3)])

axes('Position',[(0.2 + Y) (H - Z*WIN) X Z*WIN])
imagesc(Dc); axis('xy')
title('coronal');
line([0 DIM(1)],[1 1]*ORIGIN(3))
line([1 1]*ORIGIN(1),[0 DIM(3)])

axes('Position',[0.1 (H - Z*WIN - 0.1*WIN - X*WIN) Y X*WIN])
imagesc(Dt)
title 'transverse'; xlabel('y {voxels}'); ylabel('x {voxels}')
line([0 DIM(2)],[1 1]*ORIGIN(1))
line([1 1]*ORIGIN(2),[0 DIM(1)])

axes('Position',[(0.2 + Y) (H - Z*WIN - 0.1*WIN - 0.16) X 0.16])
plot(I,DS)
q      = P(max([1 (find(P == '/') + 1)]):(length(P)));
title(q,'FontSize',16,'FontWeight','Bold')
xlabel('Voxel values')
