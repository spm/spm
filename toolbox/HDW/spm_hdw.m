function spm_hdw(job)
% Warp a pair of same subject images together
%
% Very little support can be provided for the warping routine, as it
% involves optimising a very nonlinear objective function.
% Also, don't ask what the best value for the regularisation is.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_hdw.m 4148 2011-01-04 16:49:23Z guillaume $

for i=1:numel(job.data),
    run_warping(job.data(i).mov{1},job.data(i).ref{1},job.warp_opts,job.bias_opts);
end;
%_______________________________________________________________________

%_______________________________________________________________________
function run_warping(PF,PG,warp_opts,bias_opts)

VG  = spm_vol(PG);
VF  = spm_vol(PF);
reg = warp_opts.reg;
nit = warp_opts.nits;

% Load the image pair as 8 bit.
%-----------------------------------------------------------------------
bo            = bias_opts;
[VG.uint8,sf] = loaduint8(VG); % Template
VF.uint8      = bias_correction(VF,VG,bo.nits,bo.fwhm,bo.reg,bo.lmreg,sf);

% Try loading pre-existing deformation fields.  Otherwise, create
% deformation fields from uniform affine transformations.
%-----------------------------------------------------------------------
[pth,nme,ext,num] = spm_fileparts(VF.fname);
ofname = fullfile(pth,['y_' nme '.nii']);
Def    = cell(3,1);
if exist(ofname)==2,
    P      = [repmat(ofname,3,1), [',1,1';',1,2';',1,3']];
    VT     = spm_vol(P);
    if any(abs(VT(1).mat\VG.mat - eye(4))>0.00001),
        error('Incompatible affine transformation matrices.');
    end;
    Def{1} = spm_load_float(VT(1));
    Def{2} = spm_load_float(VT(2));
    Def{3} = spm_load_float(VT(3));
    spm_affdef(Def{:},inv(VF.mat));
else,
    fprintf('Generating uniform affine transformation field\n');
    Def{1} = single(1:VG.dim(1))';
    Def{1} = Def{1}(:,ones(VG.dim(2),1),ones(VG.dim(3),1));
    Def{2} = single(1:VG.dim(2));
    Def{2} = Def{2}(ones(VG.dim(1),1),:,ones(VG.dim(3),1));
    Def{3} = reshape(single(1:VG.dim(3)),1,1,VG.dim(3));
    Def{3} = Def{3}(ones(VG.dim(1),1),ones(VG.dim(2),1),:);
    spm_affdef(Def{:},VF.mat\VG.mat);
end;

% Voxel sizes
%-----------------------------------------------------------------------
vxg = sqrt(sum(VG.mat(1:3,1:3).^2))';if det(VG.mat(1:3,1:3))<0, vxg(1) = -vxg(1); end;
vxf = sqrt(sum(VF.mat(1:3,1:3).^2))';if det(VF.mat(1:3,1:3))<0, vxf(1) = -vxf(1); end;

% Do warping
%-----------------------------------------------------------------------
fprintf('Warping (iterations=%d regularisation=%g)\n', nit, reg);
spm_warp(VG.uint8,VF.uint8,Def{:},[vxg vxf],[nit,reg,1,0]);

% Convert mapping from voxels to mm
%-----------------------------------------------------------------------
spm_affdef(Def{:},VF.mat);

% Write the deformations
%-----------------------------------------------------------------------
save_def(Def,VG.mat,ofname)

return;
%_______________________________________________________________________

%_______________________________________________________________________
function [udat,sf] = loaduint8(V)
% Load data from file indicated by V into an array of unsigned bytes.

spm_progress_bar('Init',V.dim(3),...
                ['Computing max/min of ' spm_str_manip(V.fname,'t')],...
                'Planes complete');
mx = -Inf;
for p=1:V.dim(3),
    img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
    mx  = max([max(img(:)) mx]);
    spm_progress_bar('Set',p);
end;
spm_progress_bar('Init',V.dim(3),...
        ['Loading ' spm_str_manip(V.fname,'t')],...
        'Planes loaded');
sf = 255/mx;

udat = uint8(0);
udat(V.dim(1),V.dim(2),V.dim(3))=0;
for p=1:V.dim(3),
        img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
        udat(:,:,p) = round(max(min(img*sf,255),0));
        spm_progress_bar('Set',p);
end;
spm_progress_bar('Clear');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function t = transf(B1,B2,B3,T)
d2 = [size(T) 1];
t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
t  = B1*t1*B2';
return;
%_______________________________________________________________________

%_______________________________________________________________________
function udat = bias_correction(VF,VG,nits,fwhm,reg1,reg2,sf)
% This function is intended for doing bias correction prior
% to high-dimensional intra-subject registration.  Begin by
% coregistering, then run the bias correction, and then do
% the high-dimensional warping.  A version of the first image
% is returned out, that has the same bias as that of the second
% image.  This allows warping of the corrected first image, to
% match the uncorrected second image.
%
% Ideally, it would all be combined into the same model, but
% that would require quite a lot of extra work.  I should really
% sort out a good convergence criterion, and make it a proper
% Levenberg-Marquardt optimisation.

vx      = sqrt(sum(VF.mat(1:3,1:3).^2));
d       = VF.dim(1:3);
sd      = vx(1)*VF.dim(1)/fwhm; d3(1) = ceil(sd*2); krn_x   = exp(-(0:(d3(1)-1)).^2/sd.^2)/sqrt(vx(1));
sd      = vx(2)*VF.dim(2)/fwhm; d3(2) = ceil(sd*2); krn_y   = exp(-(0:(d3(2)-1)).^2/sd.^2)/sqrt(vx(2));
sd      = vx(3)*VF.dim(3)/fwhm; d3(3) = ceil(sd*2); krn_z   = exp(-(0:(d3(3)-1)).^2/sd.^2)/sqrt(vx(3));
Cbias   = kron(krn_z,kron(krn_y,krn_x)).^(-2)*prod(d)*reg1;
Cbias   = sparse(1:length(Cbias),1:length(Cbias),Cbias,length(Cbias),length(Cbias));
B3bias  = spm_dctmtx(d(3),d3(3));
B2bias  = spm_dctmtx(d(2),d3(2));
B1bias  = spm_dctmtx(d(1),d3(1));
lmRb    = speye(size(Cbias))*prod(d)*reg2;
Tbias   = zeros(d3);

ll = Inf;
spm_plot_convergence('Init','Bias Correction','- Log-likelihood','Iteration');
for subit=1:nits,

    % Compute objective function and its 1st and second derivatives
    Alpha = zeros(prod(d3),prod(d3)); % Second derivatives
    Beta  = zeros(prod(d3),1); % First derivatives
    oll   = ll;
    ll    = 0.5*Tbias(:)'*Cbias*Tbias(:);

    for z=1:VF.dim(3),
        M1  = spm_matrix([0 0 z]);
        M2  = VG.mat\VF.mat*M1;
        f1o = spm_slice_vol(VF,M1,VF.dim(1:2),0);
        f2o = spm_slice_vol(VG,M2,VF.dim(1:2),0);
        msk = (f1o==0) & (f2o==0);
        f1o(msk) = 0;
        f2o(msk) = 0;
        ro       = transf(B1bias,B2bias,B3bias(z,:),Tbias);
        msk      = abs(ro)>0.01; % false(d(1:2));

        % Use the form based on an integral for bias that is
        % far from uniform.
        f1  = f1o(msk);
        f2  = f2o(msk);
        r   = ro(msk);
        e   = exp(r);
        t1  = (f2.*e-f1);
        t2  = (f1./e-f2);
        ll  = ll + 1/4*sum(sum((t1.^2-t2.^2)./r));
        wt1 = zeros(size(f1o));
        wt2 = zeros(size(f1o));
        wt1(msk) = (2*(t1.*f2.*e+t2.*f1./e)./r + (t2.^2-t1.^2)./r.^2)/4;
        wt2(msk) = ((f2.^2.*e.^2-f1.^2./e.^2+t1.*f2.*e-t2.*f1./e)./r/2 ...
                 - (t1.*f2.*e+t2.*f1./e)./r.^2 + (t1.^2-t2.^2)./r.^3/2);

        % Use the simple symmetric form for bias close to uniform
        f1  = f1o(~msk);
        f2  = f2o(~msk);
        r   = ro(~msk);
        e   = exp(r);
        t1  = (f2.*e-f1);
        t2  = (f1./e-f2);
        ll  = ll + (sum(t1.^2)+sum(t2.^2))/4;
        wt1(~msk) = (t1.*f2.*e-t2.*f1./e)/2;
        wt2(~msk) = ((f2.*e).^2+t1.*f2.*e + (f1./e).^2+t2.*f1./e)/2;

        b3    = B3bias(z,:)';
        Beta  = Beta  + kron(b3,spm_krutil(wt1,B1bias,B2bias,0));
        Alpha = Alpha + kron(b3*b3',spm_krutil(wt2,B1bias,B2bias,1));
    end;
    spm_plot_convergence('Set',ll/prod(d));


    if subit > 1 && ll>oll,
        % Hasn't improved, so go back to previous solution
        Tbias = oTbias;
        ll    = oll;
        lmRb  = lmRb*10;
    else
        % Accept new solution
        oTbias = Tbias;
        Tbias  = Tbias(:);
        Tbias  = Tbias - (Alpha + Cbias + lmRb)\(Beta + Cbias*Tbias);
        Tbias  = reshape(Tbias,d3);
    end;
end;

udat = uint8(0);
udat(VF.dim(1),VF.dim(2),VF.dim(3)) = 0;
for z=1:VF.dim(3),
    M1 = spm_matrix([0 0 z]);
    f1 = spm_slice_vol(VF,M1,VF.dim(1:2),0);
    r  = transf(B1bias,B2bias,B3bias(z,:),Tbias);
    f1 = f1./exp(r);
    udat(:,:,z) = round(max(min(f1*sf,255),0));
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function save_def(Def,mat,fname)
% Save a deformation field as an image

dim   = [size(Def{1},1) size(Def{1},2) size(Def{1},3) 1 3];
dtype = 'FLOAT32';
off   = 0;
scale = 1;
inter = 0;
dat   = file_array(fname,dim,dtype,off,scale,inter);

N      = nifti;
N.dat  = dat;
N.mat  = mat;
N.mat0 = mat;
N.mat_intent  = 'Aligned';
N.mat0_intent = 'Aligned';
N.intent.code = 'VECTOR';
N.intent.name = 'Mapping';
N.descrip     = 'Deformation field';
create(N);
N.dat(:,:,:,1,1) = Def{1};
N.dat(:,:,:,1,2) = Def{2};
N.dat(:,:,:,1,3) = Def{3};

% Write Jacobian determinants
[pth,nme,ext,num] = spm_fileparts(fname);
jfname = fullfile(pth,['j' nme '.nii'])
dt     = spm_def2det(Def{:},N.mat);
dat    = file_array(jfname,dim(1:3),dtype,off,scale,inter);
N      = nifti;
N.dat  = dat;
N.mat  = mat;
N.mat0 = mat;
N.mat_intent  = 'Aligned';
N.mat0_intent = 'Aligned';
N.intent.code = 'NONE';
N.intent.name = 'Det';
N.descrip     = 'Jacobian Determinant';
create(N);
N.dat(:,:,:) = dt;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function theory
% This function is never called.  It simply contains my derivation of the
% objective function and its first and second derivatives for the LM
% optimisation.  Note that this would simply be summed over all voxels.
% An estimate of the variance should maybe have been included, but I'll
% leave that for a future version.

% These will be used to simplify expressions
r  = x1*p1+x2*p2;
e  = exp(r);
t1 = (f2*e-f1)
t2 = (f1/e-f2)

% This is the actual objective function
maple int((exp(-(x1*p1+x2*p2)*(1-t))*f1-exp((x1*p1+x2*p2)*t)*f2)^2/2,t=0..1)
f = 1/4*(t1^2-t2^2)/r

% These are the first derivatives
maple diff(1/4*((f1-f2*exp(x1*p1+x2*p2))^2-(f2-f1/exp(x1*p1+x2*p2))^2)/(x1*p1+x2*p2),x1)
d1 = p1*1/4*(2*(t1*f2*e+t2*f1/e)/r + (t2^2-t1^2)/r^2)

% These are the second derivatives
maple diff(diff(1/4*((f1-f2*exp(x1*p1+x2*p2))^2-(f2-f1/exp(x1*p1+x2*p2))^2)/(x1*p1+x2*p2),x1),x2)
d2 = p1*p2*(1/2*(f2^2*e^2-f1^2/e^2+t1*f2*e-t2*f1/e)/r - (t1*f2*e+t2*f1/e)/r^2 + 1/2*(t1^2-t2^2)/r^3)


% However, the above formulation has problems in limiting case of zero bias,
% so this objective function is used for these regions
f = (t1^2+t2^2)/4

% First derivs
maple diff((f1-exp(x1*p1+x2*p2)*f2)^2/4 + (f1/exp(x1*p1+x2*p2)-f2)^2/4,x1)
d1 = (p1*t1*f2*e - p1*t2*f1/e)/2

% Second derivatives
maple diff(diff((f1-exp(x1*p1+x2*p2)*f2)^2/4 + (f1/exp(x1*p1+x2*p2)-f2)^2/4,x1),x2)
d2 = p1*p2*(f2^2*e^2+t1*f2*e + f1^2/e^2+t2*f1/e)/2



