function spm_mri_sn(P,M,bb,Vox,H,I)
% Spatial normalization for T1 and T2{*} MRI data
% FORMAT spm_mri_sn(P,M,bb,Vox,H,I);
% P     - row matrix of filenames (the first should be the MRI)
% M     - filenames of the {white & gray matter smooth segmented) templates 
% bb    - bounding box in standard space {mm}
% Vox   - voxel size in standard space {mm}
% H     - H(1) - fix pitch {to specified estimate}
%         H(2) - fix roll  {to specified estimate}
%         H(3) - fix yaw   {to specified estimate}
% I     - I(1) - pitch {specified estimate}
%         I(2) - roll  {specified estimate}
%         I(3) - yaw   {specified estimate}
% 
%____________________________________________________________________________
%
% spm_mri_sn spatially normalizes structural MRI images by normalizing
% position, size and shape using a 12 parameter affine transformation.  
%
% The primary use of this stage is to:
%
% i)   Facilitate some conventional means of reporting data (e.g.
% in Talairach coordinates using a reference image that conforms to the
% Talairach space).
%
% ii)  Normalize other [PET or fMRI functional] images using the a high
% resolution MRI image.
%
% iii) Generate normalized [fuzzy] segmented MRI data for voxel-based
% morphometrics. gray*.img and white*.img.
%
% Spatial normalization is implemented by spm_mri_sn in 3 steps
%
% 1)  solving for intensity and affine transformations that 
%     match the MRI and model images.
% 2)  Applying [the inverse of] the spatial normalizing transformation 
%     to the MRI image (filenames are prefixed with a 'n')
% 3)  Applying [the inverse of] the spatial normalizing transformation
%     and the segmenting intensity transformation to the MRI image to
%     give a normalized segmented MRI (prefixed with 'seg')
%
% The image and voxel sizes of the MRI image can be different from the
% other images.  The latter would normally be PET, SPECT or fMRI data
% that are either co-registered with the MRI image or have some known spatial
% relationship to the MRI image.  This relationship is specified using the
% ORIGIN field in the header, such that the voxels specified by ORIGIN,
% in both images, are homologous and that the images have the same orientation.
% It is therefore very important that the ORIGIN field is correctly
% specified in the headers, even if they denote the same voxel (i.e. they
% are exactly congruent).
%
% The segmented images are prefixed by 'gray and white'.  The segmentation
% function is estimated simultaneously along with the spatial transformations
%
% The spatial normalization that is effected matches the object T1 or T2{*} 
% weighted MRI image and two reference, model or template images. 
% This match is in a least squares sense and the thresholding and
% spatial transformations are estimated by linearizing the problem and
% using a standard least squares solutions.
%
% The model or template images are clearly crucial and define the space into
% which the images specified are normalized.  These templates are usually
% mriWs.img and mriGs.img in the SPM directory.  They can be changed if 
% required.  The reference images should be a white and a gray matter segmented 
% image convolved with a Gaussian kernel of 8mm FWHM and should be:
%
% 1) be more 'complete' than the object image (i.e. a greater field of view)
% 2) have a complete and correct header file (including ORIGIN - usually
% taken to be the anterior commissure)
%
%
% for a complete description of this approach see Friston et al (1995)
% ref: Friston et al (1994) The spatial registration and normalization of 
% images.  HBM 2:00-00
%
%__________________________________________________________________________
% %W% %E%



%----------------------------------------------------------------------------
global CWD
cd(CWD)

% setup temporary filenames and get image parameters
%----------------------------------------------------------------------------
temp1  = [CWD '/temp1.img'];				% image string
temp2  = [CWD '/temp2.img'];				% image string
Sbox   = diff(bb) + 1;					% size of bb
Dim    = floor(diff(bb)./Vox  + 1); 			% Image size {voxels}
Vox    = Vox;						% Voxel size {mm}
Origin = -bb(1,:)./Vox + 1;				% Origin {voxels}
FWHM   = sqrt(4^2 + 8^2);				% smoothing (mm)

d      = size(P,1);
[dim vox scale bit offset origin] = spm_hread(M(1,:));	% get model parameters
[DIM VOX SCALE BIT OFFSET ORIGIN] = spm_hread(P(1,:));	% get MRI parameters
[UIM UOX SCALE BIT OFFSET URIGIN] = spm_hread(P(d,:));	% get other parameters


% open files and write headers for all specified images
%----------------------------------------------------------------------------
U      = zeros(1, size(P,1));
V      = zeros(12,size(P,1));
for i  = 1:size(P,1)
	q      = P(i,:);
	q      = q(q ~= ' ');
	V(:,i) = spm_map(q);	
	[UIM UOX SCALE TYPE OFFSET URIGIN] = spm_hread(q);
	BIT(i) = TYPE;
	d      = max([find(q == '/') 0]);
	q      = [q(1:d) 'n' q((d + 1):length(q))];
	spm_hwrite(q,Dim,Vox,SCALE,TYPE,0,Origin,['spm - normalized']);
end

% open file for segmented image
%----------------------------------------------------------------------------
q      = P(1,:);
q      = q(q ~= ' ');
d      = max([find(q == '/') 0]);
d      = [q(1:d) 'gray' q((d + 1):length(q))];
Ug     = fopen(d,'w');
spm_hwrite(q,Dim,Vox,1/255,2,0,Origin,['spm - segmented']);
d      = [q(1:d) 'white' q((d + 1):length(q))];
Uw     = fopen(d,'w');
spm_hwrite(q,Dim,Vox,1/255,2,0,Origin,['spm - segmented']);


% memory map the MRI image and model image (rCBF or gray matter-like)
%----------------------------------------------------------------------------
Vu     = spm_map(P(1,:));				% object image
Vw     = spm_map(M(1,:));				% white reference
Vg     = spm_map(M(2,:));				% gray  reference
Vs     = Vg;


% bounding boxes and centres
%----------------------------------------------------------------------------
bp      = spm_bb(P(1,:));				% MRI - bounding box
bp(1,3) = max([(bp(2,3) - 128/VOX(3)) bp(1,3)]);	% top down
bp      = round(bp);					% dimensions
bs      = [1 1 1;dim];					% model - bounding box
cents   = mean(bs);
centp   = mean(bp);

% compute matrices {dQ = 12 parameter affine transformation}
%----------------------------------------------------------------------------
Tq    = spm_matrix([0 0 0 0 0 0 UOX]);			% anisotropy of voxels
Tp    = spm_matrix([0 0 0 0 0 0 VOX]);			% anisotropy of voxels
Ts    = spm_matrix([0 0 0 0 0 0 vox]);			% anisotropy of voxels
Tq    = spm_matrix(-URIGIN.*UOX)*Tq;			% to origin
Tp    = spm_matrix(-ORIGIN.*VOX)*Tp;			% to origin
Tq    = spm_matrix([-centp.*VOX + ORIGIN.*VOX])*Tq;	% centre of bb
Tp    = spm_matrix([-centp.*VOX + ORIGIN.*VOX])*Tp;	% centre of bb
Ts    = spm_matrix([-cents.*vox ])*Ts;		% centre of bb

To    = spm_matrix((origin - mean(bs)).*vox);		% center to origin

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


% fix options
%---------------------------------------------------------------------------
if H(1); dQ(find(dQ(:,4)),:) = []; end
if H(2); dQ(find(dQ(:,5)),:) = []; end
if H(3); dQ(find(dQ(:,6)),:) = []; end


% fix scale if the field of view is limited {< 32 mm}
%---------------------------------------------------------------------------
if (DIM(1)*VOX(1) < 32); dQ(find(dQ(:,7)),:) = []; end
if (DIM(2)*VOX(2) < 32); dQ(find(dQ(:,8)),:) = []; end
if (DIM(3)*VOX(3) < 32); dQ(find(dQ(:,9)),:) = []; end

% positional vectors (sMRI) covering the object's bounding box every 4mm
%---------------------------------------------------------------------------
dx    = 4;
x     = bp(1,1):dx/VOX(1):bp(2,1); 
y     = bp(1,2):dx/VOX(2):bp(2,2);
z     = bp(1,3):dx/VOX(3):bp(2,3);

xp    = x'*ones(size(y));
yp    = ones(size(x'))*y;
zp    = zeros(size(yp));
N     = length(x)*length(y);
d     = ones(1,N);
Xp    = zeros(4, length(x)*length(y)*length(z));
sp    = [length(x) length(y) length(z)];
for i = 1:length(z)
	j  = (i - 1)*N + [1:N];
	Xp(:,j) = [xp(:)'; yp(:)'; (z(i) + zp(:)'); d];
end


% estimate segmentation function by fitting two Gaussians to histogram of 
% of voxel intensities
%---------------------------------------------------------------------------
MRI     = spm_sample_vol(Vu,Xp(1,:)',Xp(2,:)',Xp(3,:)',0);
[FS BS] = hist(MRI(MRI > mean(MRI)/8),64);
FS      = FS(:).*spm_hanning(length(FS));
FG      = spm_fit_Gaussians(BS(:),FS);

MAX1    = FG(1);
MAX2    = FG(3);
SIG1    = FG(2)/MAX1;
SIG2    = FG(4)/MAX2;


% plot histogram of voxel values and fitted Gaussians
%----------------------------------------------------------------------------
figure(3); spm_clf
r       = 0.5;
g1      = exp(-(BS - FG(1)).^2/(2*FG(2)^2));
g1      = r*g1/sum(g1);
g2      = exp(-(BS - FG(3)).^2/(2*FG(4)^2));
g2      = (1 - r)*g2/sum(g2);


subplot(2,1,1)
plot(BS,FS,BS,g1*sum(FS),'-.',BS,g2*sum(FS),'--');
title(['Fitted Gaussians - ' P(1,(P(1,:) ~= ' '))],'FontSize',16);
xlabel('Voxel intensity - segmentation functions');
ylabel('frequency')
grid on; axis square


% compute Y - segement and smooth
%---------------------------------------------------------------------------
u       = [1:size(Xp,2)];
v       = u + size(Xp,2);
Y       = zeros(2*size(Xp,2),4);

%---------------------------------------------------------------------------
sMRI    = exp(-(MRI/MAX1 - 1).^2/(2*SIG1^2));
fid     = fopen(temp1,'w');
fwrite(fid,sMRI,'float'); fclose(fid);
spm_hwrite(temp1,sp,[dx dx dx],1,16,0);
spm_smooth(temp1,temp2,FWHM);
fid     = fopen(temp2,'r');
d       = fread(fid,size(sMRI),'float'); fclose(fid);
Y(u,1)  = d;
Y(v,3)  = d;

%---------------------------------------------------------------------------
sMRI    = exp(-(MRI/MAX2 - 1).^2/(2*SIG2^2));
fid     = fopen(temp1,'w');
fwrite(fid,sMRI,'float'); fclose(fid);
spm_hwrite(temp1,sp,[dx dx dx],1,16,0);
spm_smooth(temp1,temp2,FWHM);
fid     = fopen(temp2,'r');
d       = fread(fid,size(sMRI),'float'); fclose(fid);
Y(u,2)  = d;
Y(v,4)  = d;

% initial estimate for parameters
%----------------------------------------------------------------------------
if sp(3) < 16; d = 12; else; d = 0; end
Q       = inv(spm_matrix([0 0 -d I*pi/180]));


% solve for spatial and intensity transformation parameters
% i.e. X = [-dXdQ Y].q     =>     X + dXdQ.q1 = Y.q2
%----------------------------------------------------------------------------
rCBF  = zeros(2*size(Xp,2),1);
dXdQ  = zeros(2*size(Xp,2),size(dQ,1));
GC    = zeros(12,1);
for k = 1:12
	Xq	   = inv(Ts)*Q*Tp*Xp;
	rCBF(u)    = spm_sample_vol(Vw,Xq(1,:)',Xq(2,:)',Xq(3,:)',1);
	rCBF(v)    = spm_sample_vol(Vg,Xq(1,:)',Xq(2,:)',Xq(3,:)',1);

	% compute dXdQ {implementing dQ by resampling}
	%--------------------------------------------------------------------
	for i = 1:size(dQ,1)
	 Xq        = inv(Ts)*spm_matrix(dQ(i,:))*Q*Tp*Xp;
	 dXdQ(u,i) = spm_sample_vol(Vw,Xq(1,:)',Xq(2,:)',Xq(3,:)',1) - rCBF(u);
	 dXdQ(v,i) = spm_sample_vol(Vg,Xq(1,:)',Xq(2,:)',Xq(3,:)',1) - rCBF(v);
	end

	% solve for transformation coeficients
	%-------------------------------------------------------------------
	q          = [-dXdQ Y]\rCBF;

	% update transformation matrix
	%-------------------------------------------------------------------
	for j      = 1:size(dQ,1);
		Q  = real(expm(logm(spm_matrix(dQ(j,:)))*q(j)))*Q; end
	GC(k,1)    = det(Q);
	GC(k,2)    = sqrt(sum((Q(1:3,4)*dx).^2));
end

% graph convergence
%----------------------------------------------------------------------------
subplot(2,1,2)
plot(GC)
title('convergence','FontSize',16)
xlabel('Iterations')
ylabel('translation {mm} and |Q|')
grid on


spm_print

%----------------------------------------------------------------------------
Q     = inv(Q);
Cg    = q([1 2] + size(dXdQ,2) + 2);	% select gray  matter coefficients
Cw    = q([1 2] + size(dXdQ,2) + 0);	% select white matter coefficients


% Display results of affine and quadratic transformation
%============================================================================

%----------------------------------------------------------------------------
figure(3); spm_clf
set(3,'Units','pixels')
r         = get(gcf,'Position');
r         = r(3)/r(4);
Y         = 0.36;
X         = Y*Sbox(1)/Sbox(2);
Z         = Y*Sbox(3)/Sbox(2);

% sagittal
%----------------------------------------------------------------------------
y         = ones(Dim(3),1)*[bb(1,2):Vox(2):bb(2,2)];
z         = [bb(1,3):Vox(3):bb(2,3)]'*ones(1,Dim(2));
x         = zeros(size(y));
Xt        = [x(:)'; y(:)'; z(:)'; ones(1,Dim(2)*Dim(3))];
Xs        = inv(Ts)*To*Xt;
Xp        = inv(Tp)*Q*To*Xt;

sMRI      = spm_sample_vol(Vu,Xp(1,:)',Xp(2,:)',Xp(3,:)',3);
sMRI      = reshape(sMRI,Dim(3),Dim(2));
rCBF      = spm_sample_vol(Vs,Xs(1,:)',Xs(2,:)',Xs(3,:)',3);
rCBF      = reshape(rCBF,Dim(3),Dim(2));

axes('Position',[0.1 (0.9 - Z*r) Y Z*r])
imagesc(sMRI); axis off; axis xy; title('sagittal')
line([0 Dim(2)],-bb(1,3)*[1 1]/Vox(3)); line(-[1 1]*bb(1,2)/Vox(2),[0 Dim(3)])

axes('Position',[0.1 (0.88 - 2*Z*r ) Y Z*r])
imagesc(rCBF); axis off; axis xy 
line([0 Dim(2)],-bb(1,3)*[1 1]/Vox(3)); line(-[1 1]*bb(1,2)/Vox(2),[0 Dim(3)])


% coronal
%----------------------------------------------------------------------------
x         = ones(Dim(3),1)*[bb(1,1):Vox(1):bb(2,1)];
z         = [bb(1,3):Vox(3):bb(2,3)]'*ones(1,Dim(1));
y         = zeros(size(x));
Xt        = [x(:)'; y(:)'; z(:)'; ones(1,Dim(1)*Dim(3))];
Xs        = inv(Ts)*To*Xt;
Xp        = inv(Tp)*Q*To*Xt;

sMRI      = spm_sample_vol(Vu,Xp(1,:)',Xp(2,:)',Xp(3,:)',3);
sMRI      = reshape(sMRI,Dim(3),Dim(1));
rCBF      = spm_sample_vol(Vs,Xs(1,:)',Xs(2,:)',Xs(3,:)',3);
rCBF      = reshape(rCBF,Dim(3),Dim(1));

axes('Position',[0.5 (0.9 - Z*r) X Z*r])
imagesc(sMRI); axis off; axis xy; title('coronal')
line([0 Dim(1)],-bb(1,3)*[1 1]/Vox(3)); line(Dim(1)/2*[1 1],[0 Dim(2)])

axes('Position',[0.5 (0.88 - 2*Z*r) X Z*r])
imagesc(rCBF); axis off; axis xy 
line([0 Dim(1)],-bb(1,3)*[1 1]/Vox(3)); line(Dim(1)/2*[1 1],[0 Dim(2)])


% transverse
%----------------------------------------------------------------------------
x         = [bb(1,1):Vox(1):bb(2,1)]'*ones(1,Dim(2));
y         = ones(Dim(1),1)*[bb(1,2):Vox(2):bb(2,2)];
z         = ones(size(x))*0;
Xt        = [x(:)'; y(:)'; z(:)'; ones(1,Dim(1)*Dim(2))];
Xs        = inv(Ts)*To*Xt;
Xp        = inv(Tp)*Q*To*Xt;

sMRI      = spm_sample_vol(Vu,Xp(1,:)',Xp(2,:)',Xp(3,:)',3);
sMRI      = reshape(sMRI,Dim(1),Dim(2));
rCBF      = spm_sample_vol(Vs,Xs(1,:)',Xs(2,:)',Xs(3,:)',3);
rCBF      = reshape(rCBF,Dim(1),Dim(2));

axes('Position',[0.1 (0.84 - 2*Z*r - X*r) Y X*r])
imagesc(spm_grid(sMRI)); axis off; axis xy; title 'transverse'

axes('Position',[0.1 (0.82 - 2*Z*r - 2*X*r) Y X*r])
imagesc(spm_grid(rCBF)); axis off; axis xy


% Numerical results
%----------------------------------------------------------------------------
axes('Position',[0.5 0.12 0.4 0.34],'Visible','off')
text(0,1, 'Spatial transformation','FontSize',16,'FontWeight','Bold');
text(0,0.9, P(1,:),'FontSize',10);
text(0,0.8, ['(U,V,W} = f(x,y,z)'],'FontWeight','Bold');
text(0,0.7, 'Linear {affine} matrix','FontWeight','Bold');
text(0,0.6, sprintf('U = %0.2fx + %0.2fy + %0.2fz + %0.2f',Q(1,:)));
text(0,0.5, sprintf('V = %0.2fx + %0.2fy + %0.2fz + %0.2f',Q(2,:)));
text(0,0.4, sprintf('W = %0.2fx + %0.2fy + %0.2fz + %0.2f',Q(3,:)));
set(gca,'Ylim',[0 1])


spm_print


% write images
%----------------------------------------------------------------------------
set(2,'Name','writing images');

fig   = 0;						% initialize subplot
x     = [bb(1,1):Vox(1):bb(2,1)]'*ones(1,Dim(2));	% locations {x}
y     = ones(Dim(1),1)*[bb(1,2):Vox(2):bb(2,2)];	% locations {y}
z     = ones(size(x));					% locations {z}

d     = [0:0.01:2];
SEG1  = exp(-(d - 1).^2/(2*SIG1^2));
SEG2  = exp(-(d - 1).^2/(2*SIG2^2));
MAXG  = max((Cg(1)*SEG1 + Cg(2)*SEG2));
MAXW  = max((Cw(1)*SEG1 + Cw(2)*SEG2));

for i = bb(1,3):Vox(3):bb(2,3)				% z range of bb

	Xt     = [x(:)'; y(:)'; i*z(:)'; ones(1,Dim(1)*Dim(2))];
	Xs     = inv(Ts)*To*Xt;
	Xp     = inv(Tp)*Q*To*Xt;
	Xq     = inv(Tq)*Q*To*Xt;

	% fuzzy segmentation
	%--------------------------------------------------------------------
	MRI    = spm_sample_vol(Vu,Xp(1,:)',Xp(2,:)',Xp(3,:)',1);
	MRI    = reshape(MRI,Dim(1),Dim(2));
 	rCBF   = spm_sample_vol(Vs,Xs(1,:)',Xs(2,:)',Xs(3,:)',1);
	rCBF   = reshape(rCBF,Dim(1),Dim(2));
	MASK   = (rCBF > max(rCBF(:))*0.4);
	SEG1   = exp(-(MRI/MAX1 - 1).^2/(2*SIG1^2));
	SEG2   = exp(-(MRI/MAX2 - 1).^2/(2*SIG2^2));
	SEGG   = (Cg(1)*SEG1 + Cg(2)*SEG2);
	SEGG   = SEGG.*(SEGG > 0).*MASK;
	SEGW   = (Cw(1)*SEG1 + Cw(2)*SEG2);
	SEGW   = SEGW.*(SEGW > 0);

	% write segmented images
	%--------------------------------------------------------------------
	fwrite(Ug,255*SEGG/MAXG,'uint8');
	fwrite(Uw,255*SEGW/MAXW,'uint8');

	% display images for a couple of planes
	%--------------------------------------------------------------------
 	if ((i >= 4) & (i < (4 + Vox(3)))) | (i >= 20) & (i < (20 + Vox(3)))
		if ~fig;
			figure(3); spm_clf; fig = 4; else
			fig = 7;
		end

		subplot(3,3,(fig + 0));imagesc(rot90(spm_grid(MRI)));
		axis image; axis off;
		title(['MRI @ ' num2str(i) ' mm'],'Color',[1 0 0])
		subplot(3,3,(fig + 1));imagesc(rot90(spm_grid(SEGG)));
		axis image; axis off;
		title('Fuzzy segmentation')
		subplot(3,3,(fig + 2));imagesc(rot90(spm_grid(rCBF)));
		axis image; axis off; title('template {gray}'); drawnow
	end
 
	% write the normalized MRI images
	if (i == bb(1,3)) open_mode = 'w'; else open_mode = 'a'; end
	for j = 1:size(P,1)
		q      = P(j,:);
		q      = q(q ~= ' ');
		d      = max([find(q == '/') 0]);
		q      = [q(1:d) 'n' q((d + 1):length(q))];
		if (j == 1)
		   MRI = spm_sample_vol(V(:,j),Xp(1,:)',Xp(2,:)',Xp(3,:)',1);
		else
		   MRI = spm_sample_vol(V(:,j),Xq(1,:)',Xq(2,:)',Xq(3,:)',1);
		end
		U      = fopen(q,open_mode);
		fwrite(U,MRI/V(7,j),spm_type(BIT(j)));
		fclose(U);
	end
end


% list files that have been spatially normalized
%----------------------------------------------------------------------------
subplot(3,1,1); axis off
str = '12 parameter affine (linear) normalization';
title(str,'FontSize',16,'FontWeight','Bold')

x     =  0.1;
y     =  0.9;
for i = 1:min([size(P,1) 8])
    h = text(x,y,[sprintf('%-4.0f',i) P(i,:)]);
    y = y - 0.1;
    if (i == 1); set(h,'Color',[1 0 0]); end
end

spm_print

% close, unmap and delete
%----------------------------------------------------------------------------
for i = 1:size(V,2); spm_unmap(V(:,i)); end
spm_unmap(Vu)
spm_unmap(Vs)

delete temp1.img
delete temp2.img
delete temp1.hdr
delete temp2.hdr

fclose('all');
