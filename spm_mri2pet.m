function spm_mri2pet
% co-registration of (T1-weighted) MRI to PET
% FORMAT spm_mri2pet;
% 
%____________________________________________________________________________
%
% spm_mri2pet is a rough mri to pet co-registration that perfoms a fuzzy
% segmentation of the mri image, convolves it and then matches it to the
% pet image in the usual way (using least squares and first order Taylor
% expansions).
% spm_mri2pet will prompt for the requisite filenames and write the
% registered mri image *.img to r*.img
%
% This routine accepts T1-weighted mri images that have not been scalp 
% edited or pre-processed in any way.
%
% The pet data should be the average of any time-series and should be
% devoid of blanks or missing data.
%
% for a description of this approach see Friston et al (1994)
% ref: Friston et al (1994) The spatial registration and normalization of 
% images. HBM 0:00-00
%
%__________________________________________________________________________
% %W% %E%


%----------------------------------------------------------------------------
global CWD
cd(CWD)

P      = spm_get(1,'.img','Select PET image');
M      = spm_get(1,'.img','Select MRI image');

%----------------------------------------------------------------------------
figure(2); clf; set(2,'Name','Thankyou','Pointer','Watch'); drawnow

%----------------------------------------------------------------------------
temp1  = [CWD '/temp1.img'];				% image string
temp2  = [CWD '/temp2.img'];				% image string
temp3  = [CWD '/temp3.img'];				% image string
temp4  = [CWD '/temp4.img'];				% image string
FWHM   = 8;						% smoothing (mm)
Vox    = [2 2 2];


% Open files and write headers
%----------------------------------------------------------------------------
Vu     = spm_map(P);
Vv     = spm_map(M);	
q      = max([find(M == '/') 0]);
q      = [M(1:q) 'r' M((q + 1):length(M))];
U      = fopen(q(q ~= ' '),'w');
[DIM VOX SCALE TYPE OFFSET ORIGIN] = spm_hread(P);
[dim vox scale bits offset origin] = spm_hread(M);
spm_hwrite(q,DIM,VOX,scale,bits,0,ORIGIN,['mri - co-registered']);



% segment, smooth and memory map
%============================================================================

% fuzzy segment MRI image
%----------------------------------------------------------------------------
bs      = spm_bb(M);
bs(1,:) = bs(1,:) - 16/vox(1);
bs(2,:) = bs(2,:) + 16/vox(1);
ds      = round(diff(bs).*vox/2);			% size in 2mm voxels

To      = spm_matrix([-bs(1,:) 0 0 0 1 1 1]);		% corner of bb
To      = spm_matrix([0 0 0 0 0 0 vox/2])*To;		% anisotropy of voxels

[FS BS] = spm_hist(M(M ~= ' '),dim,bits);
BS      = BS*scale;
FS      = FS(:).*hanning(length(FS));
FG      = spm_fit_Gaussians(BS(:),FS);

MAX1    = FG(1);
MAX2    = FG(3);
SIG1    = FG(2)/MAX1;
SIG2    = FG(4)/MAX2;
MAX     = min([MAX1 MAX2]);
if MAX == MAX1; SIG = SIG1/MAX1; else; SIG = SIG2/MAX2; end

% plot histogram of voxel values and fitted Gaussians
%----------------------------------------------------------------------------
figure(3); spm_clf
r       = 0.5;
g1      = exp(-(BS - FG(1)).^2/(2*FG(2)^2));
g1      = r*g1/sum(g1);
g2      = exp(-(BS - FG(3)).^2/(2*FG(4)^2));
g2      = (1 - r)*g2/sum(g2);

plot(BS,FS,BS,g1*sum(FS),'-.',BS,g2*sum(FS),'--');
title(['Fitted Gaussians - ' M(M ~= ' ')],'FontSize',16);
xlabel('Voxel intensity - segmentation functions');
ylabel('frequency')
grid on; axis square

spm_print

% write segmented image
%----------------------------------------------------------------------------
Q       = fopen(temp1,'w');
for i   = 1:ds(3)
	d = spm_slice_vol(Vv,inv(spm_matrix([0 0 -i])*To),ds(1:2),0);
	fwrite(Q,(255*exp(-(d/MAX - 1).^2/(2*(1/8)^2))),'uint8');
end
fclose(Q);


% fuzzy segment PET image
%----------------------------------------------------------------------------
bp      = spm_bb(P);
bp      = bp + [-8/VOX(1) -8/VOX(1) 0;8/VOX(1) 8/VOX(1) 0];
dp      = round(diff(bp));				% size in 2mm voxels
Tn      = spm_matrix([-bp(1,:) 0 0 0 1 1 1]);		% corner of bb
d       = spm_slice_vol(Vu,inv(spm_matrix([0 0 -dp(3)/2])*Tn),dp(1:2),0);
q       = max(d(:));
Q       = fopen(temp2,'w');
for i   = 1:dp(3)
	d = spm_slice_vol(Vu,inv(spm_matrix([0 0 -i])*Tn),dp(1:2),0);
	fwrite(Q,(255*(1 + tanh(d/q*6 - 3))/2),'uint8');
end
fclose(Q);

spm_hwrite(temp1,ds,Vox,1,2,0);
spm_hwrite(temp2,dp,VOX,1,2,0);
spm_hwrite(temp3,ds,Vox,1,2,0);
spm_hwrite(temp4,dp,VOX,1,2,0);

spm_smooth(temp1,temp3,sqrt(FWHM.^2 + 8^2));		% smooth mri
spm_smooth(temp2,temp4,FWHM);				% smooth pet

Vs    = spm_map(temp3);					% mri image
Vp    = spm_map(temp4);					% pet image


% center, bounding box and margins
%----------------------------------------------------------------------------
bp    = [1 1 1;dp];					% pet - bounding box
bs    = [1 1 1;ds];					% mri - bounding box

% compute matrices {dQ = 12 affine transformations, dW = 9 second order terms}
%----------------------------------------------------------------------------
Tp    = spm_matrix([-mean(bp) 0 0 0 1 1 1]);		% centre of bb
Ts    = spm_matrix([-mean(bs) 0 0 0 1 1 1]);		% centre of bb
Tp    = spm_matrix([0 0 0 0 0 0 Vp(4:6)'])*Tp;		% anisotropy of voxels
Ts    = spm_matrix([0 0 0 0 0 0 Vs(4:6)'])*Ts;		% anisotropy of voxels


r     = 0.1*pi/180;					% dQ (rotation) radians
t     = 0.1;						% dQ (translation) mm
s     = 1.01;						% dQ (scale) 
a     = 0.01;						% dQ (affine) 
dQ    = [t 0 0 0 0 0 1 1 1 0 0 0;			% unit transformations
 	 0 t 0 0 0 0 1 1 1 0 0 0;
	 0 0 t 0 0 0 1 1 1 0 0 0;
         0 0 0 r 0 0 1 1 1 0 0 0;
         0 0 0 0 r 0 1 1 1 0 0 0;
         0 0 0 0 0 r 1 1 1 0 0 0];


% positional vectors (PET) covering the object's bounding box
%---------------------------------------------------------------------------

x     = bp(1,1):4/Vp(4):bp(2,1); 
y     = bp(1,2):4/Vp(5):bp(2,2);
z     = bp(1,3):6/Vp(6):bp(2,3);

xp    = x'*ones(size(y));
yp    = ones(size(x'))*y;
zp    = zeros(size(yp));
N     = length(x)*length(y);
d     = ones(1,N);
Xp    = zeros(4, length(x)*length(y)*length(z));
for i = 1:length(z)
	j  = (i - 1)*N + [1:N];
	Xp(:,j) = [xp(:)'; yp(:)'; (z(i) + zp(:)'); d];
end


% compute PET {X} the object image and effects of movement (dXdQ and dXdW)
%---------------------------------------------------------------------------
PET   = zeros(size(Xp,2),1);
SPM   = zeros(size(Xp,2),1);
dXdQ  = zeros(size(Xp,2),size(dQ,1));


% compute PET the object image
%---------------------------------------------------------------------------
PET   = spm_sample_vol(Vp,Xp(1,:)',Xp(2,:)',Xp(3,:)',1);

% compute dXdQ {implementing dQ by resampling}
%---------------------------------------------------------------------------
for i = 1:size(dQ,1)
	Xq        = inv(Tp)*spm_matrix(dQ(i,:))*Tp*Xp;
	dXdQ(:,i) = spm_sample_vol(Vp,Xq(1,:)',Xq(2,:)',Xq(3,:)',1) - PET;
end

% initial estimate for parameters
%----------------------------------------------------------------------------
Q     = spm_matrix([0 0 -16 16*pi/180 0 0]);
G     = (PET + mean(PET));

for k = 1:16

	Xq	  = inv(Ts)*inv(Q)*Tp*Xp;
	SPM       = spm_sample_vol(Vs,Xq(1,:)',Xq(2,:)',Xq(3,:)',1);
	SPM       = SPM.*G;
	SPM       = [SPM.^0 SPM.^1 SPM.^2];

	% solve for transformation coeficients
	%-------------------------------------------------------------------
	q         = [-dXdQ SPM]\PET;

	% update transformation matrices
	%-------------------------------------------------------------------
	for j     = 1:size(dQ,1);
		Q = real(expm(logm(spm_matrix(dQ(j,:)))*q(j)))*Q;
	end

end

% final trasnfromation matrix
%----------------------------------------------------------------------------
A      = inv(To)*inv(Ts)*inv(Q)*Tp*Tn;				% MRI -> PET

% Display results of affine and quadratic transformation
%============================================================================
figure(3); spm_clf
load Split; colormap(split)
set(3,'Units','pixels')

bb     = (bp - [1; 1]*mean(bp)).*([1; 1]*VOX);		% bb in mm
sb     = diff(bb) + 1;					% size of bb
Dim    = floor(diff(bb)./Vox  + 1); 			% dimensions (i j and k}

r      = get(gcf,'Position');
r      = r(3)/r(4);
Y      = 0.34;
X      = Y*sb(1)/sb(2);
Z      = Y*sb(3)/sb(2);

% sagittal
%----------------------------------------------------------------------------
y         = ones(Dim(3),1)*[bb(1,2):Vox(2):bb(2,2)];
z         = [bb(1,3):Vox(3):bb(2,3)]'*ones(1,Dim(2));
x         = zeros(size(y));
Xt        = [x(:)'; y(:)'; z(:)'; ones(1,Dim(2)*Dim(3))];
Xp        = inv(Tp)*Xt;
Xi        = inv(Ts)*inv(Q)*Xt;
Xs        = inv(To)*inv(Ts)*inv(Q)*Xt;

PET       = spm_sample_vol(Vp,Xp(1,:)',Xp(2,:)',Xp(3,:)',1);
PET       = reshape(PET,Dim(3),Dim(2));
SPM       = spm_sample_vol(Vs,Xi(1,:)',Xi(2,:)',Xi(3,:)',1);
SPM       = reshape(SPM,Dim(3),Dim(2));
MRI       = spm_sample_vol(Vv,Xs(1,:)',Xs(2,:)',Xs(3,:)',1);
MRI       = reshape(MRI,Dim(3),Dim(2));

axes('Position',[0.1 (0.9 - Z*r) Y Z*r])
image(PET*64/max(PET(:))); axis off; axis xy; title('sagittal')
line([0 Dim(2)],-bb(1,3)*[1 1]/Vox(3)); line(-[1 1]*bb(1,2)/Vox(2),[0 Dim(3)])

axes('Position',[0.1 (0.88 - 2*Z*r ) Y Z*r])
image(SPM*64/max(SPM(:))); axis off; axis xy 
line([0 Dim(2)],-bb(1,3)*[1 1]/Vox(3)); line(-[1 1]*bb(1,2)/Vox(2),[0 Dim(3)])

axes('Position',[0.1 (0.86 - 3*Z*r ) Y Z*r])
image(MRI*64/max(MRI(:))); axis off; axis xy 
line([0 Dim(2)],-bb(1,3)*[1 1]/Vox(3)); line(-[1 1]*bb(1,2)/Vox(2),[0 Dim(3)])



% coronal
%----------------------------------------------------------------------------
x         = ones(Dim(3),1)*[bb(1,1):Vox(1):bb(2,1)];
z         = [bb(1,3):Vox(3):bb(2,3)]'*ones(1,Dim(1));
y         = zeros(size(x));
Xt        = [x(:)'; y(:)'; z(:)'; ones(1,Dim(1)*Dim(3))];
Xp        = inv(Tp)*Xt;
Xi        = inv(Ts)*inv(Q)*Xt;
Xs        = inv(To)*inv(Ts)*inv(Q)*Xt;

PET       = spm_sample_vol(Vp,Xp(1,:)',Xp(2,:)',Xp(3,:)',1);
PET       = reshape(PET,Dim(3),Dim(1));
SPM       = spm_sample_vol(Vs,Xi(1,:)',Xi(2,:)',Xi(3,:)',1);
SPM       = reshape(SPM,Dim(3),Dim(1));
MRI       = spm_sample_vol(Vv,Xs(1,:)',Xs(2,:)',Xs(3,:)',1);
MRI       = reshape(MRI,Dim(3),Dim(1));

axes('Position',[0.5 (0.9 - Z*r) X Z*r])
image(PET*64/max(PET(:))); axis off; axis xy; title('coronal')
line([0 Dim(1)],-bb(1,3)*[1 1]/Vox(3)); line(Dim(1)/2*[1 1],[0 Dim(2)])

axes('Position',[0.5 (0.88 - 2*Z*r) X Z*r])
image(SPM*64/max(SPM(:))); axis off; axis xy 
line([0 Dim(1)],-bb(1,3)*[1 1]/Vox(3)); line(Dim(1)/2*[1 1],[0 Dim(2)])

axes('Position',[0.5 (0.86 - 3*Z*r) X Z*r])
image(MRI*64/max(MRI(:))); axis off; axis xy 
line([0 Dim(1)],-bb(1,3)*[1 1]/Vox(3)); line(Dim(1)/2*[1 1],[0 Dim(2)])


% transverse
%----------------------------------------------------------------------------
x         = [bb(1,1):Vox(1):bb(2,1)]'*ones(1,Dim(2));
y         = ones(Dim(1),1)*[bb(1,2):Vox(2):bb(2,2)];
z         = ones(size(x)) - 8;
Xt        = [x(:)'; y(:)'; z(:)'; ones(1,Dim(1)*Dim(2))];
Xp        = inv(Tn)*inv(Tp)*Xt;
Xs        = inv(To)*inv(Ts)*inv(Q)*Xt;

PET       = spm_sample_vol(Vu,Xp(1,:)',Xp(2,:)',Xp(3,:)',1);
PET       = reshape(PET,Dim(1),Dim(2));
SPM       = spm_sample_vol(Vv,Xs(1,:)',Xs(2,:)',Xs(3,:)',1);
SPM       = reshape(SPM,Dim(1),Dim(2));


axes('Position',[0.1 (0.82 - 3*Z*r - X*r) Y X*r])
image(PET*64/max(PET(:)));
axis off;  title('transverse')

axes('Position',[0.1 (0.80 - 3*Z*r - 2*X*r) Y X*r])
image(SPM*64/max(SPM(:)));
axis off;

% Numerical results
%----------------------------------------------------------------------------
axes('Position',[0.5 0.06 0.4 0.34],'Visible','off')
text(0,1, 'MRI -> PET registration','FontSize',16,'FontWeight','Bold');
text(0,0.9, ['{MRI} ' M],'FontSize',10);
text(0,0.8, ['{PET} ' P],'FontSize',10);
text(0,0.7, '(U V W) = f (x y z)','FontWeight','Bold');
text(0,0.6, sprintf('U = %0.2fx + %0.2fy + %0.2fz + %0.2f',A(1,:)));
text(0,0.5, sprintf('V = %0.2fx + %0.2fy + %0.2fz + %0.2f',A(2,:)));
text(0,0.4, sprintf('W = %0.2fx + %0.2fy + %0.2fz + %0.2f',A(3,:)));
text(0,0.3, 'Transformation matrix','FontWeight','Bold');
set(gca,'Ylim',[0 1])

spm_print

% display merged data sets
%----------------------------------------------------------------------------
figure(3); spm_clf

bb     = (2*bs - [1; 1]*mean(2*bs)).*([1; 1]*vox);		% bb in mm
sb     = diff(bb) + 1;					% size of bb
Dim    = floor(diff(bb)./Vox  + 1); 			% dimensions (i j and k}

r      = get(gcf,'Position');
r      = r(3)/r(4);
Y      = 0.38;
X      = Y*sb(1)/sb(2);
Z      = Y*sb(3)/sb(2);

% sagittal
%----------------------------------------------------------------------------
y         = ones(Dim(3),1)*[bb(1,2):Vox(2):bb(2,2)];
z         = [bb(1,3):Vox(3):bb(2,3)]'*ones(1,Dim(2));
x         = zeros(size(y));
Xt        = [x(:)'; y(:)'; z(:)'; ones(1,Dim(2)*Dim(3))];
Xp        = inv(Tp)*Xt;
Xs        = inv(To)*inv(Ts)*inv(Q)*Xt;

PET       = spm_sample_vol(Vp,Xp(1,:)',Xp(2,:)',Xp(3,:)',1);
PET       = reshape(PET,Dim(3),Dim(2));
MRI       = spm_sample_vol(Vv,Xs(1,:)',Xs(2,:)',Xs(3,:)',1);
MRI       = reshape(MRI,Dim(3),Dim(2));

axes('Position',[0.1 (0.9 - Z*r) Y Z*r])
image(spm_merge(128*PET/max(PET(:)),64*MRI/max(MRI(:))));
axis off; axis xy; title('sagittal')


% coronal
%----------------------------------------------------------------------------
x         = ones(Dim(3),1)*[bb(1,1):Vox(1):bb(2,1)];
z         = [bb(1,3):Vox(3):bb(2,3)]'*ones(1,Dim(1));
y         = zeros(size(x));
Xt        = [x(:)'; y(:)'; z(:)'; ones(1,Dim(1)*Dim(3))];
Xp        = inv(Tp)*Xt;
Xs        = inv(To)*inv(Ts)*inv(Q)*Xt;

PET       = spm_sample_vol(Vp,Xp(1,:)',Xp(2,:)',Xp(3,:)',1);
PET       = reshape(PET,Dim(3),Dim(1));
MRI       = spm_sample_vol(Vv,Xs(1,:)',Xs(2,:)',Xs(3,:)',1);
MRI       = reshape(MRI,Dim(3),Dim(1));

axes('Position',[0.52 (0.9 - Z*r) X Z*r])
image(spm_merge(128*PET/max(PET(:)),64*MRI/max(MRI(:))));
axis off; axis xy; title('coronal')

% transverse
%----------------------------------------------------------------------------
x         = [bb(1,1):Vox(1):bb(2,1)]'*ones(1,Dim(2));
y         = ones(Dim(1),1)*[bb(1,2):Vox(2):bb(2,2)];
z         = ones(size(x)) - 8;
Xt        = [x(:)'; y(:)'; z(:)'; ones(1,Dim(1)*Dim(2))];
Xp        = inv(Tn)*inv(Tp)*Xt;
Xs        = inv(To)*inv(Ts)*inv(Q)*Xt;

PET       = spm_sample_vol(Vu,Xp(1,:)',Xp(2,:)',Xp(3,:)',1);
PET       = reshape(PET,Dim(1),Dim(2));
MRI       = spm_sample_vol(Vv,Xs(1,:)',Xs(2,:)',Xs(3,:)',1);
MRI       = reshape(MRI,Dim(1),Dim(2));

axes('Position',[0.1 (0.80 - Z*r - X*r) Y X*r])
image(spm_merge(128*PET/max(PET(:)),64*MRI/max(MRI(:))));
axis off; title('transverse')


% text
%----------------------------------------------------------------------------
axes('Position',[0.52 0.16 0.5 0.34],'Visible','off')
text(0,1, 'MRI -> PET registration','FontSize',16,'FontWeight','Bold');
text(0,0.9, ['{MRI} ' M],'FontSize',10);
text(0,0.8, ['{PET} ' P],'FontSize',10);
set(gca,'Ylim',[0 1])



spm_print


% write image transformed mri image
%============================================================================

%----------------------------------------------------------------------------
x      = [1:DIM(1)]'*ones(1,DIM(2));
y      = ones(DIM(1),1)*[1:DIM(2)];
z      = ones(size(x));

for i = 1:DIM(3)

	Xq    = A*[x(:)'; y(:)'; i*z(:)'; ones(1,DIM(1)*DIM(2))];

	% write transformed images
	%--------------------------------------------------------------------
	MRI    = spm_sample_vol(Vv,Xq(1,:)',Xq(2,:)',Xq(3,:)',3);
	MRI    = reshape(MRI,DIM(1),DIM(2));
	MRI    = MRI.*(MRI > 0);
	fwrite(U,MRI/scale,spm_type(bits));
end


% close, unmap and delete
%----------------------------------------------------------------------------
spm_unmap_vol(Vu)
spm_unmap_vol(Vv)
spm_unmap_vol(Vp)
spm_unmap_vol(Vs)

delete temp1.img
delete temp1.hdr
delete temp2.img
delete temp2.hdr
delete temp3.img
delete temp3.hdr
delete temp4.img
delete temp4.hdr

fclose('all');
set(2,'Name',' ','Pointer','Arrow');

