
% rendering of regional effects [SPM{Z}] on orthogonal sections
% FORMAT spm_sections
%____________________________________________________________________________
%
% spm_sections is called by spm_results and uses variables in working
% memory to create three orthogonal sections though a background image.
% Regional foci from the selected SPM{Z} are rendered on this image.
%
% Although the SPM{Z} adopts the neurological convention (left = left)
% the rendered images follow the same convention as the original data.
%
%__________________________________________________________________________
% %W% Karl Friston %E%



% get image on which to render
%----------------------------------------------------------------------------
spms   = spm_get(1,'.img','select an image for rendering');
set([Finter,Fgraph],'Pointer','Watch');

% memory map background image and create transformation matrix {A}
%----------------------------------------------------------------------------
L      = spm_XYZreg('RoundCoords',L,V);
Vs     = spm_vol(spms);					% memory mapped
Nx     = round(V(1)*V(4));				% SPM size (mm)
Ny     = round(V(2)*V(5));				% SPM size (mm)
O      = round(V(3)*V(6));				% SPM size (mm)
J      = round(V(7)*V(4));				% corner of SPM (mm)
R      = round(V(8)*V(5));				% corner of SPM (mm)
I      = round(V(9)*V(6));				% corner of SPM (mm)
A      = spm_matrix([J R 0])*Vs.mat;			% re-center to corner


% extract data from SPM{t} [sagittal {Ts} coronal {Tc} transverse {Tt}]
%----------------------------------------------------------------------------
Q      = find(abs(XYZ(1,:) - L(1)) < V(4)/2);
Ts     = sparse((XYZ(3,Q) + I)/V(6),(XYZ(2,Q) + R)/V(5),t(Q),O/V(6),Ny/V(5));
Ts     = spm_resize(full(Ts),O,Ny);

Q      = find(abs(XYZ(2,:) - L(2)) < V(5));
Tc     = sparse((XYZ(3,Q) + I)/V(6),(XYZ(1,Q) + J)/V(4),t(Q),O/V(6),Nx/V(4));
Tc     = spm_resize(full(Tc),O,Nx);

Q      = find(abs(XYZ(3,:) - L(3)) < V(6));
Tt     = sparse((XYZ(1,Q) + J)/V(4),(XYZ(2,Q) + R)/V(5),t(Q),Nx/V(4),Ny/V(5));
Tt     = spm_resize(full(Tt),Nx,Ny);


% get background slices and combine
%----------------------------------------------------------------------------
D      = [0 0 1 I;0 1 0 0;-1 0 0 (J + L(1));0 0 0 1]*A;
D      = spm_slice_vol(Vs,inv(D),[O Ny],1);
Ds     = D/max(D(:));

D      = [0 0 1 I;1 0 0 0;0 -1 0 (R + L(2));0 0 0 1]*A;
D      = spm_slice_vol(Vs,inv(D),[O Nx],1);
Dc     = D/max(D(:));

D      = [1 0 0 0;0 1 0 0;0 0 1 -L(3);0 0 0 1]*A;
D      = spm_slice_vol(Vs,inv(D),[Nx Ny],1);
Dt     = D/max(D(:));


%-Delete previous axis and their pagination controls (if any)
%-----------------------------------------------------------------------
spm_results_ui('ClearPane',Fgraph,'RNP');


% configure {128 level} colormap
%----------------------------------------------------------------------------
load Split
colormap(split)

d      = max([Ts(:); Tc(:); Tt(:)]);
D      = length(colormap)/2;
Q      = Ts(:) > U; Ts = Ts(Q)/d; Ds(Q) = 1 + Ts; Ts = D*Ds;
Q      = Tc(:) > U; Tc = Tc(Q)/d; Dc(Q) = 1 + Tc; Tc = D*Dc;
Q      = Tt(:) > U; Tt = Tt(Q)/d; Dt(Q) = 1 + Tt; Tt = D*Dt;


% compute axes to correct for anisotropy of voxels and (normalized) window
%----------------------------------------------------------------------------
set(Fgraph,'Units','pixels')
WIN    = get(gcf,'Position');
WIN    = WIN(3)/WIN(4);
Y      = 0.36;
X      = Y*Nx/Ny;
Z      = Y*O/Ny;

% render activation foci on background image
%----------------------------------------------------------------------------
axes('Position',[0.1 (0.46 - Z*WIN) Y Z*WIN])
image(Ts)
axis image; axis('xy'); axis off; 	title 'sagittal'
line([0 Ny],([I I] + L(3)))
line(([R R] + L(2)),[0 O])

axes('Position',[(0.2 + Y) (0.46 - Z*WIN) X Z*WIN])
image(Tc)
axis image; axis('xy'); axis off; 	title 'coronal'
line([0 Nx],([I I] + L(3)))
line(([J J] + L(1)),[0 O])

axes('Position',[0.1 (0.46 - Z*WIN - 0.1*WIN - X*WIN) Y X*WIN])
image(Tt)
axis image; axis off;
title 'transverse'
line([0 Ny],([J J] + L(1)))
line(([R R] + L(2)),[0 Nx])


% colorbar
%----------------------------------------------------------------------------
q     = [(0.2 + Y) (0.46 - Z*WIN - 0.1*WIN - X*WIN) 0.02 X*WIN];
axes('Position',q,'Visible','off')
image([0 d/32],[0 d],[1:D]' + D)

if SPMZ str ='Z value'; end;
if SPMF str = 'F-value'; end;

axis xy; ylabel(str);
set(gca,'XTickLabel',[])


% reset pointer (and x locations is necessary)
%----------------------------------------------------------------------------
set([Finter,Fgraph],'Pointer','Arrow')
