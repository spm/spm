function spm_sn(P,qmean,spms,bb,Vox,H,I)
% Spatial normalization {PET -> model or template}
% FORMAT spm_sn(P,q,Q,bb,Vox,H,[I]);
% P     - row matrix of filenames of aligned images
% q     - vector of indices indicating files to be used in the average
% Q     - filename of template image
% bb    - bounding box for standard space {mm}
% Vox   - voxel size for standard space {mm}
% H     - vector of options:
% 	H(1) - fix pitch to specified estimate
% 	H(2) - 6 parameter 3D-nonlinear (quadratic)
% 	H(3) - 2D-nonlinear (Fourier-like basis functions) for each slice
%       H(4) - number of parameters of the affine component
%               6  - parameters: - A rigid body transformation
%               7  - parameters: - A size and position only transformation
%               9  - parameters: - A size and orthogonal scaling transformation
%               12 - parameters: - A full affine transformation
%
% I     - initial pitch estimate {degrees}
% 
%____________________________________________________________________________
%
% Spatial normalization is implemented by spm_sn in 4 steps
%
% 1)  taking the mean of several co-registered images
% 2)  solving for an affine/[quadratic] transformation, matching mean to a 
%     model
% 3)  slice-based nonlinear normalization (using Fourier-like basis functions)
% 4)  writing the normalized image (filenames are prefixed with a 'n')
%
% The spatial normalization that is effected matches the mean object image
% (specified in terms of the filenames to be averaged) and a reference, model
% or template image.  This match is in a least squares sense and the spatial
% nonlinear transformations are estimated by linearizing the problem and
% using a standard explicit least squares solution.
%
% There are two levels of nonlinear normalization that can be used.  The first
% is large scale and fully 3-dimensional (the Quadratic component).  The second
% effects transformations at a finer spatial scale and works slice by slice.
% The slice-based transformations include an intensity transformation that
% modulates the template in an optimal way.  This allows for subjects with
% large infarcts and the like.  The 'corrected' template is displayed 
% for a couple of slices.
%
% The model or template image is clearly crucial and defines the space into
% which the images specified are normalized.  This template is usually the 
% image called spm.img in the SPM directory {SWD}.  It can be changed to any 
% reference required.  The reference image should:
%
% 1) be more 'complete' than the object images (i.e. a greater field of view)
% 2) be of the same modality and resolution as the object image
% 3) be relatively smooth
% 4) have a complete and correct header file (including ORIGIN)
%
% The smoothness requirement relates to the first order approximations used
% in the linearization of the problem.
%
% for a complete description of this approach see Friston et al (1994)
% ref: Friston et al (1994) The spatial registration and normalization of 
% images.  HBM 0:00-00
%
%__________________________________________________________________________
% %W% %E%

%----------------------------------------------------------------------------
global PRINTSTR GRID SWD CWD
cd(CWD)

%----------------------------------------------------------------------------
temp1  = [CWD '/temp1.img'];				% mean  image string
temp2  = [CWD '/temp2.img'];				% mean  image string
sb     = diff(bb) + 1;					% size of bb
Dim    = floor(diff(bb)./Vox  + 1); 			% dimensions (i j and k}
origin = -bb(1,:)./Vox + 1;				% origin {voxels}
FWHM   = 8;						% smoothing (mm)


% Open files and write headers
%----------------------------------------------------------------------------
U      = zeros(1, size(P,1));
V      = zeros(12,size(P,1));
for i  = 1:size(P,1)
	q      = P(i,:);
	q      = q(q ~= ' ');
	V(:,i) = spm_map(q);	
	[DIM VOX SCALE TYPE OFFSET] = spm_hread(q);
	TYP(i) = TYPE;
	d      = max([find(q == '/') 0]);
	q      = [q(1:d) 'n' q((d + 1):length(q))];
	U(i)   = fopen(q,'w');
	spm_hwrite(q,Dim,Vox,SCALE,TYPE,0,origin,['spm - normalized']);
end

% check for consistency of image size and voxel size
%----------------------------------------------------------------------------
if ~(all(all(~diff(V([1:6],:)')))) & (size(P,1) > 1)
	error('images are not compatible'); end

if ~(all(V([4:6],1)))
	error('voxels have zero volume! '); end


% average selected images, smooth and memory map
%----------------------------------------------------------------------------
spm_mean(prod(DIM),TYP(1),temp1,P(qmean,:));		% average image
spm_hwrite(temp1,DIM,VOX,1,TYP(1),OFFSET);		% write header
spm_smooth(temp1,temp2,FWHM);				% smooth average
Vu    = spm_map(temp1);					% mean object image
Vp    = spm_map(temp2);					% mean object (smoothed)
Vs    = spm_map(spms);					% model reference image


% center, bounding box and margins
%----------------------------------------------------------------------------
[DIM d d d d ORIGIN] = spm_hread(spms);			% get DIM & ORIGIN
bp    = spm_bb(temp1);					% bb - bounding box
bs    = [1 1 1;DIM];					% bb - bounding box


% compute matrices {dQ = 12 affine transformations, dW = 9 second order terms}
%----------------------------------------------------------------------------
Tp    = spm_matrix([-mean(bp) 0 0 0 1 1 1]);		% centre of bb
Ts    = spm_matrix([-mean(bs) 0 0 0 1 1 1]);		% centre of bb
Tp    = spm_matrix([0 0 0 0 0 0 Vp(4:6)'])*Tp;		% anisotropy of voxels
Ts    = spm_matrix([0 0 0 0 0 0 Vs(4:6)'])*Ts;		% anisotropy of voxels
To    = spm_matrix((ORIGIN - mean(bs)).*Vs(4:6)');	% center to ORIGIN

r     = 0.1*pi/180;					% dQ (rotation) radians
t     = 0.1;						% dQ (translation) mm
s     = 1.01;						% dQ (scale) 
a     = 0.01;						% dQ (affine) 
dQ    = [t 0 0 0 0 0 1 1 1 0 0 0;			% unit transformations
 	 0 t 0 0 0 0 1 1 1 0 0 0;
	 0 0 t 0 0 0 1 1 1 0 0 0;
         0 0 0 r 0 0 1 1 1 0 0 0;
         0 0 0 0 r 0 1 1 1 0 0 0;
         0 0 0 0 0 r 1 1 1 0 0 0;
         0 0 0 0 0 0 s 1 1 0 0 0;
         0 0 0 0 0 0 1 s 1 0 0 0;
         0 0 0 0 0 0 1 1 s 0 0 0;
         0 0 0 0 0 0 1 1 1 a 0 0;
         0 0 0 0 0 0 1 1 1 0 a 0;
         0 0 0 0 0 0 1 1 1 0 0 a];

dW    = [1 0 0 0 0 0 0 0 0;				% quadratic terms
	 0 1 0 0 0 0 0 0 0;				% omitting W = f(Z^2)
	 0 0 1 0 0 0 0 0 0;
	 0 0 0 1 0 0 0 0 0;
	 0 0 0 0 1 0 0 0 0;
	 0 0 0 0 0 1 0 0 0];

dW    = 0.001*dW;

% number of paramters for affine component
%---------------------------------------------------------------------------
dQ    = dQ([1:H(4)],:);
if H(4) == 7; dQ(7,:) = [0 0 0 0 0 0 s s s 0 0 0]; end

% fix pitch option
%---------------------------------------------------------------------------
if H(1); dQ(find(dQ(:,4)),:) = []; end

% quadratic option
%---------------------------------------------------------------------------
if ~H(2); dW = []; end



% positional vectors (PET) covering the object's bounding box
%---------------------------------------------------------------------------
Xp    = [];
x     = bp(1,1):4/Vp(4):bp(2,1); 
y     = bp(1,2):4/Vp(5):bp(2,2);
z     = (bp(1,3) + FWHM/V(6)):8/Vp(6):(bp(2,3) - FWHM/V(6));
M     = length(x);
N     = length(y);
for i = z
	xp = x'*ones(size(y));
	yp = ones(size(x'))*y;
	zp = i*ones(size(yp));
	Xp = [Xp [xp(:)'; yp(:)'; zp(:)'; ones(1,M*N)]];
end


% remove redundant parts {i.e less than mean/8)
%---------------------------------------------------------------------------
PET   = spm_sample_vol(Vu,Xp(1,:)',Xp(2,:)',Xp(3,:)',1);
d     = PET > mean(PET)/8;
Xp    = Xp(:,d);

% compute PET {X} the object image and effects of movement (dXdQ and dXdW)
%---------------------------------------------------------------------------
PET   = zeros(size(Xp,2),1);
SPM   = zeros(size(Xp,2),1);
dXdQ  = zeros(size(Xp,2),size(dQ,1));
dXdW  = zeros(size(Xp,2),size(dW,1));


% compute PET the object image
%---------------------------------------------------------------------------
PET   = spm_sample_vol(Vp,Xp(1,:)',Xp(2,:)',Xp(3,:)',1);

% compute dXdQ {implementing dQ by resampling}
%---------------------------------------------------------------------------
for i = 1:size(dQ,1)
	Xq        = inv(Tp)*spm_matrix(dQ(i,:))*Tp*Xp;
	dXdQ(:,i) = spm_sample_vol(Vp,Xq(1,:)',Xq(2,:)',Xq(3,:)',1) - PET;
end

% compute dXdW {implementing dW by resampling}
%---------------------------------------------------------------------------
for i = 1:size(dW,1)
	Xq        = Tp*Xp;
	Xq(1:3,:) = Xq(1:3,:) + reshape(dW(i,:),3,3)*Xq(1:3,:).^2;
	Xq	  = inv(Tp)*Xq;
	dXdW(:,i) = spm_sample_vol(Vp,Xq(1,:)',Xq(2,:)',Xq(3,:)',1) - PET;
end


% initial estimate for parameters
%----------------------------------------------------------------------------
Q     = spm_matrix([0 0 0 I*pi/180 0 0]);
W     = zeros(1,9);


% compute SPM the reference image and solve for Q and W
% i.e. PET = [-dXdQ -dXdW SPM]*q
% NB an approximation to the inverse of the quadratic component is used here
% but it appears to be a good approximation over the ranges encountered
%----------------------------------------------------------------------------
m     = max(PET);
d     = m/4;
g     = d:d:m;
G     = zeros(size(PET,1),size(g,2));
for i = 1:length(g)
	G(:,i) = exp(-(g(i) - PET).^2/(2*d^2)); end


for k = 1:8
	Xq        = Tp*Xp;
	Xq(1:3,:) = Xq(1:3,:) - reshape(W,3,3)*Xq(1:3,:).^2;
	Xq	  = inv(Ts)*inv(Q)*Xq;
	SPM       = spm_sample_vol(Vs,Xq(1,:)',Xq(2,:)',Xq(3,:)',1);
	SPM       = SPM*ones(1,size(G,2)).*G;

	% solve for transformation coeficients
	%-------------------------------------------------------------------
	q         = [-dXdQ -dXdW SPM]\PET;

	% update transformation matrices
	%-------------------------------------------------------------------
	for j     = 1:size(dQ,1);
		Q = real(expm(logm(spm_matrix(dQ(j,:)))*q(j)))*Q; end
	for j     = 1:size(dW,1)
		W = W + dW(j,:)*q(j + size(dQ,1)); 		  end

	% impose constraints on transformations
	%-------------------------------------------------------------------
	if Q(3,3) > 1.16 | Q(3,3) < 1/1.16
		d         = find(dQ(:,9) ~=1);
		dQ(d,:)   = [];
		dXdQ(:,d) = [];
	end
	d         = [];
	for j     = 1:size(dW,1)
		if abs(W(find(dW(j,:)))) > 0.0005; d = [d j]; end
	end
	if ~isempty(d)
		dW(d,:)   = [];
		dXdW(:,d) = [];
	end
end

W      = reshape(W,3,3);

% Display results of affine and quadratic transformation
%============================================================================
figure(3); spm_clf
set(3,'Units','pixels')
r      = get(gcf,'Position');
r      = r(3)/r(4);
Y      = 0.36;
X      = Y*sb(1)/sb(2);
Z      = Y*sb(3)/sb(2);

% sagittal
%----------------------------------------------------------------------------
y         = ones(Dim(3),1)*[bb(1,2):Vox(2):bb(2,2)];
z         = [bb(1,3):Vox(3):bb(2,3)]'*ones(1,Dim(2));
x         = zeros(size(y));
Xt        = [x(:)'; y(:)'; z(:)'; ones(1,Dim(2)*Dim(3))];
Xs        = inv(Ts)*To*Xt;
Xp        = Q*To*Xt;
Xp(1:3,:) = Xp(1:3,:) + W*Xp(1:3,:).^2;
Xp        = inv(Tp)*Xp;

PET       = spm_sample_vol(Vu,Xp(1,:)',Xp(2,:)',Xp(3,:)',1);
PET       = reshape(PET,Dim(3),Dim(2));
SPM       = spm_sample_vol(Vs,Xs(1,:)',Xs(2,:)',Xs(3,:)',1);
SPM       = reshape(SPM,Dim(3),Dim(2));

axes('Position',[0.1 (0.9 - Z*r) Y Z*r])
imagesc(PET); axis off; axis xy; title('sagittal')
line([0 Dim(2)],-bb(1,3)*[1 1]/Vox(3)); line(-[1 1]*bb(1,2)/Vox(2),[0 Dim(3)])

axes('Position',[0.1 (0.88 - 2*Z*r ) Y Z*r])
imagesc(SPM); axis off; axis xy 
line([0 Dim(2)],-bb(1,3)*[1 1]/Vox(3)); line(-[1 1]*bb(1,2)/Vox(2),[0 Dim(3)])


% coronal
%----------------------------------------------------------------------------
x         = ones(Dim(3),1)*[bb(1,1):Vox(1):bb(2,1)];
z         = [bb(1,3):Vox(3):bb(2,3)]'*ones(1,Dim(1));
y         = zeros(size(x));
Xt        = [x(:)'; y(:)'; z(:)'; ones(1,Dim(1)*Dim(3))];
Xs        = inv(Ts)*To*Xt;
Xp        = Q*To*Xt;
Xp(1:3,:) = Xp(1:3,:) + W*Xp(1:3,:).^2;
Xp        = inv(Tp)*Xp;

PET       = spm_sample_vol(Vu,Xp(1,:)',Xp(2,:)',Xp(3,:)',1);
PET       = reshape(PET,Dim(3),Dim(1));
SPM       = spm_sample_vol(Vs,Xs(1,:)',Xs(2,:)',Xs(3,:)',1);
SPM       = reshape(SPM,Dim(3),Dim(1));

axes('Position',[0.5 (0.9 - Z*r) X Z*r])
imagesc(PET); axis off; axis xy; title('coronal')
line([0 Dim(1)],-bb(1,3)*[1 1]/Vox(3)); line(Dim(1)/2*[1 1],[0 Dim(2)])

axes('Position',[0.5 (0.88 - 2*Z*r) X Z*r])
imagesc(SPM); axis off; axis xy 
line([0 Dim(1)],-bb(1,3)*[1 1]/Vox(3)); line(Dim(1)/2*[1 1],[0 Dim(2)])


% transverse
%----------------------------------------------------------------------------
x         = [bb(1,1):Vox(1):bb(2,1)]'*ones(1,Dim(2));
y         = ones(Dim(1),1)*[bb(1,2):Vox(2):bb(2,2)];
z         = zeros(size(x));
Xt        = [x(:)'; y(:)'; z(:)'; ones(1,Dim(1)*Dim(2))];
Xs        = inv(Ts)*To*Xt;
Xp        = Q*To*Xt;
Xp(1:3,:) = Xp(1:3,:) + W*Xp(1:3,:).^2;
Xp        = inv(Tp)*Xp;

PET       = spm_sample_vol(Vu,Xp(1,:)',Xp(2,:)',Xp(3,:)',1);
PET       = reshape(PET,Dim(1),Dim(2));
SPM       = spm_sample_vol(Vs,Xs(1,:)',Xs(2,:)',Xs(3,:)',1);
SPM       = reshape(SPM,Dim(1),Dim(2));

axes('Position',[0.1 (0.84 - 2*Z*r - X*r) Y X*r])
imagesc(spm_grid(PET)); axis off; axis xy; title 'transverse'

axes('Position',[0.1 (0.82 - 2*Z*r - 2*X*r) Y X*r])
imagesc(spm_grid(SPM)); axis off; axis xy


% Numerical results
%----------------------------------------------------------------------------
axes('Position',[0.5 0.16 0.4 0.34],'Visible','off')
text(0,1, 'Spatial transformation','FontSize',16,'FontWeight','Bold');
text(0,0.9, [P(1,:) ' etc'],'FontSize',10);
text(0,0.8, ['(u,v,w} = f(x,y,z)'],'FontWeight','Bold');
text(0,0.7, 'Linear {affine} component','FontWeight','Bold');
text(0,0.6, sprintf('U = %0.2fx + %0.2fy + %0.2fz + %0.2f',Q(1,:)));
text(0,0.5, sprintf('V = %0.2fx + %0.2fy + %0.2fz + %0.2f',Q(2,:)));
text(0,0.4, sprintf('W = %0.2fx + %0.2fy + %0.2fz + %0.2f',Q(3,:)));
text(0,0.3, 'Nonlinear {quadratic} component','FontWeight','Bold');
s = sprintf('u = U + (%0.2fU^2 + %0.2fV^2)/1000',W(1,[1:2])*1000);
text(0,0.2,s,'FontSize',10);
s = sprintf('v = V + (%0.2fU^2 + %0.2fV^2)/1000',W(2,[1:2])*1000);
text(0,0.1,s,'FontSize',10);
s = sprintf('w = W + (%0.2fU^2 + %0.2fV^2)/1000',W(3,[1:2])*1000);
text(0,0.0,s,'FontSize',10);
set(gca,'Ylim',[0 1])

spm_print



% nonlinear transformation using Fourier-like basis functions
%============================================================================

% create basis functions BASIS and global scope
%----------------------------------------------------------------------------
M     = Dim(1);						% x - {voxels}
N     = Dim(2);						% y - {voxels}
BASIS = spm_basis(M,N);					% basis functions
k     = size(BASIS,2);					% number
q     = zeros(3*k,1);					% initialize q
fig   = 0;						% initialize subplot
x     = [bb(1,1):Vox(1):bb(2,1)]'*ones(1,Dim(2));	% locations {x}
y     = ones(Dim(1),1)*[bb(1,2):Vox(2):bb(2,2)];	% locations {y}
z     = ones(size(x));					% locations {z}


for i = bb(1,3):Vox(3):bb(2,3)				% z range of bb
	Xt    = [x(:)'; y(:)'; i*z(:)'; ones(1,Dim(1)*Dim(2))];
	Xs    = inv(Ts)*To*Xt;				% model vectors
	SPM   = spm_sample_vol(Vs,Xs(1,:)',Xs(2,:)',Xs(3,:)',1);
	SPM   = reshape(SPM,Dim(1),Dim(2));		% model image
	B     = zeros(2*k,1);				% initalize B
	if H(3)						% conditional skip
	    for j = 1:3					% iterate

		% compute objects and reference inages (PET and SPM)
		%------------------------------------------------------------
		Xp        = Q*To*Xt;
		dx        = BASIS*B([1:k]    );
		dy        = BASIS*B([1:k] + k);
		Xp(1:3,:) = Xp(1:3,:) + W*Xp(1:3,:).^2;
		Xp(1:2,:) = Xp(1:2,:) + [dx(:)'; dy(:)'];
		Xp        = inv(Tp)*Xp;
		PET       = spm_sample_vol(Vp,Xp(1,:)',Xp(2,:)',Xp(3,:)',1);
		PET       = reshape(PET,Dim(1),Dim(2));
		D         = (SPM(:) > max(SPM(:))/6) | (PET(:) > max(PET(:))/6);
		MOD       = SPM.*(~(~PET));

		% compute dXdQ 
		%------------------------------------------------------------ 
		[dy dx]   = gradient(PET);
		dXdB      = [  (-dx(D)*ones(1,k)).*BASIS(D,:),...
	   		       (-dy(D)*ones(1,k)).*BASIS(D,:),...
          		       (MOD(D)*ones(1,k)).*BASIS(D,:)];

		% least squares solution
		%-----------------------------------------------------------
	    	q         = dXdB\PET(D);
	   	q         = q.*finite(q); 
	    	B         = B + q(1:(2*k));
	    end
	end
	Xp        = Q*To*Xt;
	dx        = BASIS*B([1:k]    );
	dy        = BASIS*B([1:k] + k);
	Xp(1:3,:) = Xp(1:3,:) + W*Xp(1:3,:).^2;
	Xp(1:2,:) = Xp(1:2,:) + [dx(:)'; dy(:)'];
	Xp        = inv(Tp)*Xp;


	% display images for a couple of planes
	%--------------------------------------------------------------------
 	if ((i >= 0) & (i < Vox(3))) | (i >= 20) & (i < (20 + Vox(3)))
		if ~fig;
			figure(3); spm_clf; fig = 4; else
			fig = 7; end

		PET   = spm_sample_vol(Vp,Xp(1,:)',Xp(2,:)',Xp(3,:)',1);
		PET   = reshape(PET,Dim(1),Dim(2));
		OBJ   = spm_sample_vol(Vu,Xp(1,:)',Xp(2,:)',Xp(3,:)',1);
		OBJ   = reshape(OBJ,Dim(1),Dim(2));
		SPM   = SPM.*(~(~PET));
	   	if H(3)
			SPM(:) = SPM(:).*(BASIS*q([1:k] + 2*k)); end

		subplot(3,3,(fig + 0));imagesc(rot90(spm_grid(OBJ)));
		axis image; axis off; title 'unsmoothed'
		subplot(3,3,(fig + 1));imagesc(rot90(spm_grid(PET)));
		axis image; axis off;
		title(['nonlinear @ ' num2str(i) ' mm'],'Color',[1 0 0])
		subplot(3,3,(fig + 2));imagesc(rot90(spm_grid(SPM)));
		axis image; axis off; title 'template'; drawnow
	end
 
	% write transformed images
	%--------------------------------------------------------------------
	for j = 1:length(U)
		PET    = spm_sample_vol(V(:,j),Xp(1,:)',Xp(2,:)',Xp(3,:)',1);;
		fwrite(U(j),PET/V(7,j),spm_type(TYP(j)));
	end
end


% display files spatially normalized
%----------------------------------------------------------------------------
subplot(3,1,1); axis off
if ( H(2) &  H(3)); str = 'Nonlinear spatial normalization'; end
if ( H(2) & ~H(3)); str = 'Affine and (3D) quadratic normalization'; end
if (~H(2) &  H(3)); str = 'Affine and (2D) nonlinear normalization'; end
if (~H(2) & ~H(3))
	str = sprintf('%0.0i parameter affine normalization',H(4));
end

title(str,'FontSize',16,'FontWeight','Bold')

x     =  0.1;
y     =  0.9;
for i = 1:size(P,1) 
    h = text(x,y,[sprintf('%-4.0f',i) P(i,:)]);
    y = y - 0.1;
    if qmean(i); set(h,'Color',[1 0 0]); end
end

spm_print

% close, unmap and delete
%----------------------------------------------------------------------------
for i = 1:size(V,2); spm_unmap(V(:,i)); end
spm_unmap(Vu)
spm_unmap(Vp)
spm_unmap(Vs)

delete temp1.img
delete temp1.hdr
delete temp2.img
delete temp2.hdr

fclose('all');

