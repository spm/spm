
% Rendering of regional effects [SPM{Z/F}] on transverse sections
% FORMAT spm_transverse
%_______________________________________________________________________
%
% spm_transverse is called by spm_results and uses variable in working
% memory to create three transverse sections though a background image.
% Regional foci from the selected SPM{Z/F} are rendered on this image.
%
% Although the SPM{Z} adopts the neurological convention (left = left)
% the rendered images follow the same convention as the original data.
%
%_______________________________________________________________________
% %W% Karl Friston %E%


% get the image on which to render
%-----------------------------------------------------------------------
spms   = spm_get(1,'.img','select an image for rendering');
Finter = spm_figure('FindWin','Interactive');
Fgraph = spm_figure('FindWin','Graphics');
set([Finter,Fgraph],'Pointer','Watch');

[d d d d d origin] = spm_hread(spms);
if V(3) == 1
	origin = [0 0 0]; L(3) = V(6); end

% memory map the background image and create transformation matrix {A}
%-----------------------------------------------------------------------
L      = V(4:6)'.*round(L./V(4:6)');
if FLIP; XYZ(1,:) = -XYZ(1,:); L(1) = -L(1); end	% left = right
Vs     = spm_map(spms);					% memory mapped
M      = round(V(1)*V(4));				% SPM size (mm)
N      = round(V(2)*V(5));				% SPM size (mm)
J      = round(V(7)*V(4));				% corner of SPM (mm)
Z      = round(V(8)*V(5));				% corner of SPM (mm)
A      = spm_matrix(-origin);				% center on origin
A      = spm_matrix([0 0 0 0 0 0 Vs(4:6)'])*A;		% anisotropy
A      = spm_matrix([J Z 0])*A;				% re-center to corner



% extract data from SPM{t} [at one plane separation]
%-----------------------------------------------------------------------
Q      = find(abs(XYZ(3,:) - L(3)) < V(6));
T      = sparse((XYZ(1,Q) + J)/V(4),(XYZ(2,Q) + Z)/V(5),t(Q),M/V(4),N/V(5));
Tc     = spm_resize(full(T),M,N);

if V(3) > 1
    Q  = find(abs(XYZ(3,:) - L(3) - V(6)) < V(6));
    T  = sparse((XYZ(1,Q) + J)/V(4),(XYZ(2,Q) + Z)/V(5),t(Q),M/V(4),N/V(5));
    Ts = spm_resize(full(T),M,N);

    Q  = find(abs(XYZ(3,:) - L(3) + V(6)) < V(6));
    T  = sparse((XYZ(1,Q) + J)/V(4),(XYZ(2,Q) + Z)/V(5),t(Q),M/V(4),N/V(5));
    Tt = spm_resize(full(T),M,N);
end

% get background slices and combine
%-----------------------------------------------------------------------
D      = [1 0 0 0;0 1 0 0;0 0 1 -L(3);0 0 0 1]*A;
D      = spm_slice_vol(Vs,inv(D),[M N],1);
Dc     = D/max(D(:));

if V(3) > 1
    D  = [1 0 0 0;0 1 0 0;0 0 1 (-L(3) - V(6));0 0 0 1]*A;
    D  = spm_slice_vol(Vs,inv(D),[M N],1);
    Ds = D/max(D(:));

    D  = [1 0 0 0;0 1 0 0;0 0 1 (-L(3) + V(6));0 0 0 1]*A;
    D  = spm_slice_vol(Vs,inv(D),[M N],1);
    Dt = D/max(D(:));
end

% delete previous axis
%-----------------------------------------------------------------------
figure(Fgraph)
subplot(2,1,2); delete(gca), spm_figure('DeletePageControls')

% configure {128 level} colormap
%-----------------------------------------------------------------------
load Split
colormap(split)

d      = max([Ts(:); Tc(:); Tt(:)]);
D      = length(colormap)/2;
Q      = Tc(:) > U; Tc = Tc(Q)/d; Dc(Q) = 1 + Tc; Tc = D*Dc;

if V(3) > 1
    Q  = Ts(:) > U; Ts = Ts(Q)/d; Ds(Q) = 1 + Ts; Ts = D*Ds;
    Q  = Tt(:) > U; Tt = Tt(Q)/d; Dt(Q) = 1 + Tt; Tt = D*Dt;
end

% render activation foci on background images
%-----------------------------------------------------------------------
if V(3) > 1
	subplot(2,4,6)
	image(rot90(spm_grid(Tc)))
	axis image; axis off;
	if FLIP
		title(sprintf('z = %0.0fmm {left = right}',L(3)));
	else
		title(sprintf('z = %0.0fmm',L(3)))
	end
	line(([J J] + L(1)),[0 N])
	line([0 M],((N - Z)*[1 1] - L(2)))

	subplot(2,4,5)
	image(rot90(spm_grid(Tt)))
	axis image; axis off;
	title(sprintf('z = %0.0fmm',(L(3) - V(6))))
	line(([J J] + L(1)),[0 N])
	line([0 M],((N - Z)*[1 1] - L(2)))

	subplot(2,4,7)
	image(rot90(spm_grid(Ts)))
	axis image; axis off; title(sprintf('z = %0.0fmm',(L(3) + V(6))))
	line(([J J] + L(1)),[0 N])
	line([0 M],((N - Z)*[1 1] - L(2)))
else
	axes('position', [0.3 0.1 0.4 0.3])
	image(rot90(spm_grid(Tc)))
	axis image; axis off;
	if FLIP
		title(sprintf('z = %0.0fmm {left = right}',L(3)));
	else
		title(sprintf('z = %0.0fmm',L(3)))
	end
	line(([J J] + L(1)),[0 N])
	line([0 M],((N - Z)*[1 1] - L(2)))
end

% colorbar
%-----------------------------------------------------------------------
u     = get(gca,'Position');
axes('position', [(u(1) + u(3) + 0.1) u(2) 0.01 u(3)])
image([0 d/32],[0 d],[1:D]' + D)
if SPMZ str = 'Z value'; end;
if SPMF str = 'F-value'; end;

axis xy; ylabel(str);

set(gca,'XTickLabels',[])


% unmap and reset pointer (and x locations if necessary)
%-----------------------------------------------------------------------
if FLIP; XYZ(1,:) = -XYZ(1,:); L(1) = -L(1); end	% left = right
spm_unmap(Vs);
set([Finter,Fgraph],'Pointer','Arrow')
