function vbm_preproc(P)
% A quick and dirty script for VBM preprocessing
if nargin==0,
    P=spm_get(Inf,'*.IMAGE','Select T1 weighted images');
end;
for i=1:size(P,1),
    do_vbm(P(i,:));
end;
return;
%=======================================================================

%=======================================================================
function do_vbm(P)
[pth,nam,ext]=fileparts(P);
disp(['************* ' nam ' *************']);
p = spm_preproc(P);
savefields(['vbm_' nam '.mat'],p)
%p = load(['vbm_' nam '.mat']);
preproc_apply(p)
return;
%=======================================================================

%=======================================================================
function savefields(fnam,p)
if length(p)>1, error('Can''t save fields.'); end;
fn = fieldnames(p);
for i=1:length(fn),
    eval([fn{i} '= p.' fn{i} ';']);
end;
save(fnam,fn{:});
return;
%=======================================================================

%=======================================================================
function preproc_apply(p)
[pth,nam,ext]=fileparts(p.image.fname);
T    = p.Twarp;
bsol = p.Tbias;
d2   = size(T);
V    = p.image;
d    = V.dim(1:3);
B    = p.priors;
%M   = B(1).mat\p.Affine*V(1).mat;
M    = p.Affine*V(1).mat;

[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3  = 1:d(3);
d3  = [size(bsol) 1];
B1  = spm_dctmtx(d(1),d2(1));
B2  = spm_dctmtx(d(2),d2(2));
B3  = spm_dctmtx(d(3),d2(3));
D1  = spm_dctmtx(d(1),d2(1),'diff');
D2  = spm_dctmtx(d(2),d2(2),'diff');
D3  = spm_dctmtx(d(3),d2(3),'diff');
bB3 = spm_dctmtx(d(3),d3(3),x3);
bB2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
bB1 = spm_dctmtx(d(1),d3(1),x1(:,1));

mg  = p.mg;
mn  = p.mn;
vr  = p.vr;
K   = length(p.mg);
Kb  = length(p.priors)+1;

dat = uint8(0);
dat(d(1),d(2),d(3)) = 0;
dat = {dat,dat,dat};

lkp = []; for k=1:Kb, lkp = [lkp ones(1,p.ngaus(k))*k]; end;
b0  = spm_load_priors(B);

disp('Creating Segmented...');
M = B(1).mat\p.Affine*V(1).mat;
for z=1:length(x3),
    [t1,t2,t3] = defs(T,z,B1,B2,B3,x1,x2,x3,M);
    f  = spm_sample_vol(V,x1,x2,o*x3(z),0);
    msk=find(f<p.thresh);
    cr = exp(transf(bB1,bB2,bB3(z,:),bsol)).*f;
    q  = zeros([d(1:2) Kb]);
    bg = ones(d(1:2));
    bt = zeros([d(1:2) Kb]);
    for k1=1:Kb,
        bt(:,:,k1) = spm_sample_priors(b0{k1},t1,t2,t3,k1==Kb);
    end;
    b = zeros([d(1:2) K]);
    for k=1:K,
        b(:,:,k) = bt(:,:,lkp(k))*mg(k);
    end;
    s = sum(b,3);
    for k=1:K,
        p1 = exp((cr-mn(k)).^2/(-2*vr(k)))/sqrt(2*pi*vr(k)+eps);
        q(:,:,lkp(k)) = q(:,:,lkp(k)) + p1.*b(:,:,k)./s;
    end;
    sq = sum(q,3)+eps;
    for k1=1:3,
        tmp = q(:,:,k1);
        tmp(msk) = 0;
        dat{k1}(:,:,z) = uint8(round(255 * tmp./sq));
    end;
end;

%disp('Cleaning...');
%[dat{:}] = clean_gwc(dat{:});

disp('Writing Segmented...');
fn = {['rawGM_' nam '.img'],['rawWM_' nam '.img'],['rawCSF_' nam '.img']};
for k1=1:3,
    VO        = V;
    VO.fname  = fn{k1};
    VO.dim(4) = spm_type('uint8');
    VO.pinfo  = [1/255 0 0]';
    VO        = spm_create_vol(VO);
    for z=1:VO.dim(3),
        VO    = spm_write_plane(VO,double(dat{k1}(:,:,z))/255,z);
    end;
    VO        = spm_close_vol(VO);
end;


od   = p.priors(1).dim(1:3);
M1   = p.priors(1).mat;
M2   = p.image(1).mat;
y1   = single(0);
y1(d(1),d(2),d(3)) = 0;
y2   = y1;
y3   = y1;
M    = p.Affine*V(1).mat;
disp('Creating Deformations...');
for z=1:d(3),
    [y1(:,:,z),y2(:,:,z),y3(:,:,z)] = defs(T,z,B1,B2,B3,x1,x2,x3,M);
end;

disp('Inverting Deformations...');
[y1,y2,y3] = spm_invdef(y1,y2,y3,od,inv(M1),M2);
spm_affdef(y1,y2,y3,inv(M2));


ovx      = abs(det(B(1).mat(1:3,1:3)))^(1/3);
fwhm     = max(ovx./sqrt(sum(V(1).mat(1:3,1:3).^2))-1,0.1);
dat{1}   = decimate(dat{1},fwhm);
dat{2}   = decimate(dat{2},fwhm);
tmp      = B(1).mat/V(1).mat;
tmp(4,:) = [0 0 0 1];
dt       = spm_def2det(y1,y2,y3,tmp);

disp('Warping...');

dat1 = zeros(size(y1));
dat1 = {dat1,dat1,dat1};
for z=1:size(dat1{1},3),
    for k1=1:3,
        tmp = spm_sample_vol(dat{k1},double(y1(:,:,z)),double(y2(:,:,z)),double(y3(:,:,z)),1);
        dat1{k1}(:,:,z)  = 1/255*tmp.*double(dt(:,:,z));
        %dat1{k1}(:,:,z) = 1/255*tmp;
    end;
end;
fn = {['GM_' nam '.img'],['WM_' nam '.img'],['CSF_' nam '.img']};
for k1=1:3,
    VO = struct('fname',fn{k1},'dim',[size(dat1{1}) 4],'mat',B(1).mat);
    spm_write_vol(VO,dat1{k1});
end;

dat1 = zeros(size(y1));
for z=1:size(dat1,3);
    tmp = spm_sample_vol(V,double(y1(:,:,z)),double(y2(:,:,z)),double(y3(:,:,z)),1);
    dat1(:,:,z) = tmp;
end;
VO = struct('fname',['w' nam '.img'],'dim',[size(dat1) 2],'mat',B(1).mat);
spm_write_vol(VO,dat1);
return;
%=======================================================================

%=======================================================================
function [x1,y1,z1] = defs(sol,z,B1,B2,B3,x0,y0,z0,M)
x1a = x0    + transf(B1,B2,B3(z,:),sol(:,:,:,1));
y1a = y0    + transf(B1,B2,B3(z,:),sol(:,:,:,2));
z1a = z0(z) + transf(B1,B2,B3(z,:),sol(:,:,:,3));
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
return;
%=======================================================================

%=======================================================================
function d = jacobian(T,j,B1,B2,B3,D1,D2,D3)
tx = get_2Dtrans(T(:,:,:,1),B3,j);
ty = get_2Dtrans(T(:,:,:,2),B3,j);
tz = get_2Dtrans(T(:,:,:,3),B3,j);

j11 = D1*tx*B2' + 1; j12 = B1*tx*D2';     j13 = B1*get_2Dtrans(T(:,:,:,1),D3,j)*B2';
j21 = D1*ty*B2';     j22 = B1*ty*D2' + 1; j23 = B1*get_2Dtrans(T(:,:,:,2),D3,j)*B2';
j31 = D1*tz*B2';     j32 = B1*tz*D2';     j33 = B1*get_2Dtrans(T(:,:,:,3),D3,j)*B2' + 1;

% The determinant of the Jacobian reflects relative volume changes.
%------------------------------------------------------------------
d  = (j11.*(j22.*j33-j23.*j32) - j21.*(j12.*j33-j13.*j32) + j31.*(j12.*j23-j13.*j22));
%=======================================================================

%=======================================================================
function T2 = get_2Dtrans(T3,B,j)
d   = [size(T3) 1 1 1];
tmp = reshape(T3,d(1)*d(2),d(3));
T2  = reshape(tmp*B(j,:)',d(1),d(2));
return;
%=======================================================================

%=======================================================================
function t = transf(B1,B2,B3,T)
d2 = [size(T) 1];
t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
t  = B1*t1*B2';
return;
%=======================================================================

%=======================================================================
function dat = decimate(dat,fwhm)
% Convolve the volume in memory (fwhm in voxels).
lim = ceil(2*fwhm);
s  = fwhm/sqrt(8*log(2));
x  = [-lim(1):lim(1)]; x = spm_smoothkern(fwhm(1),x); x  = x/sum(x);
y  = [-lim(2):lim(2)]; y = spm_smoothkern(fwhm(2),y); y  = y/sum(y);
z  = [-lim(3):lim(3)]; z = spm_smoothkern(fwhm(3),z); z  = z/sum(z);
i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
spm_conv_vol(dat,dat,x,y,z,-[i j k]);
return;
%=======================================================================

%=======================================================================
function [g,w,c] = clean_gwc(g,w,c)
b    = w;
b(1) = w(1);

% Build a 3x3x3 seperable smoothing kernel
%-----------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

% Erosions and conditional dilations
%-----------------------------------------------------------------------
niter = 32;
spm_progress_bar('Init',niter,'Extracting Brain','Iterations completed');
for j=1:niter,
        if j>2, th=0.15; else th=0.6; end; % Dilate after two its of erosion.
        for i=1:size(b,3),
                gp = double(g(:,:,i));
                wp = double(w(:,:,i));
                bp = double(b(:,:,i))/255;
                bp = (bp>th).*(wp+gp);
                b(:,:,i) = uint8(round(bp));
        end;
        spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
        spm_progress_bar('Set',j);
end;
th = 0.05;
for i=1:size(b,3),
        gp       = double(g(:,:,i))/255;
        wp       = double(w(:,:,i))/255;
        cp       = double(c(:,:,i))/255;
        bp       = double(b(:,:,i))/255;
        bp       = ((bp>th).*(wp+gp))>th;
        g(:,:,i) = uint8(round(255*gp.*bp./(gp+wp+cp+eps)));
        w(:,:,i) = uint8(round(255*wp.*bp./(gp+wp+cp+eps)));
        c(:,:,i) = uint8(round(255*(cp.*bp./(gp+wp+cp+eps)+cp.*(1-bp))));
end;
spm_progress_bar('Clear');
return;
%=======================================================================

%=======================================================================

