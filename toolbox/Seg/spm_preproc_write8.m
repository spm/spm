function cls = spm_preproc_write(res,bf,tc)
% Write out VBM preprocessed data
% FORMAT spm_preproc_write(p,opts)
% p    - results from spm_prep2sn
% opts - writing options.  A struct containing these fields:
%        biascor - write bias corrected image
%        GM      - flags for which images should be written
%        WM      - similar to GM
%        CSF     - similar to GM
%____________________________________________________________________________
% Copyright (C) 2008 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_preproc_write8.m 1230 2008-03-19 14:23:20Z john $

tpm = res.tpm;
if ~isstruct(tpm) || ~isfield(tpm, 'bg'),
    tpm = spm_load_priors8(tpm);
end

lkp = res.lkp;
Kb  = max(lkp);

[pth,nam,ext]=fileparts(res.image(1).fname);
ind  = res.image(1).n;
d    = res.image(1).dim(1:3);

[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3  = 1:d(3);

N   = numel(res.image);
chan(N) = struct('B1',[],'B2',[],'B3',[],'T',[],'Nc',[],'Nf',[],'ind',[]);
for n=1:N,
    d3  = [size(res.Tbias{n}) 1];
    chan(n).B3 = spm_dctmtx(d(3),d3(3),x3);
    chan(n).B2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
    chan(n).B1 = spm_dctmtx(d(1),d3(1),x1(:,1));
    chan(n).T  = res.Tbias{n};

    [pth1,nam1,ext1] = fileparts(res.image(n).fname);
    chan(n).ind      = res.image(n).n;

    if bf(n,2),
        chan(n).Nc      = nifti;
        chan(n).Nc.dat  = file_array(fullfile(pth1,['m', nam1, ext1]),...
                                     res.image(n).dim(1:3),...
                                     [spm_type('float32') spm_platform('bigend')],...
                                     0,1,0);
        chan(n).Nc.mat  = res.image(n).mat;
        chan(n).Nc.mat0 = res.image(n).mat;
        chan(n).Nc.descrip = 'Bias corrected';
        create(chan(n).Nc);
    end

    if bf(n,1),
        chan(n).Nf      = nifti;
        chan(n).Nf.dat  = file_array(fullfile(pth1,['BiasField_', nam1, ext1]),...
                                     res.image(n).dim(1:3),...
                                     [spm_type('float32') spm_platform('bigend')],...
                                     0,1,0);
        chan(n).Nf.mat  = res.image(n).mat;
        chan(n).Nf.mat0 = res.image(n).mat;
        chan(n).Nf.descrip = 'Estimated Bias Field';
        create(chan(n).Nf);
    end
end

do_cls = false;
tiss(Kb) = struct('Nt',[]);
cls = cell(1,Kb);
for k1=1:Kb,
    if tc(k1,2),
        cls{k1} = zeros(d(1:3),'uint8');
        do_cls  = true;
    end
    if tc(k1,1),
        tiss(k1).Nt      = nifti;
        tiss(k1).Nt.dat  = file_array(fullfile(pth,['c', num2str(k1), nam, ext1]),...
                                      res.image(n).dim(1:3),...
                                      [spm_type('uint8') spm_platform('bigend')],...
                                      0,1/255,0);
        tiss(k1).Nt.mat  = res.image(n).mat;
        tiss(k1).Nt.mat0 = res.image(n).mat;
        tiss(k1).Nt.descrip = ['Tissue class ' num2str(k1)];
        create(tiss(k1).Nt);
        do_cls = true;
    end;
end

prm  = [3 3 3 0 0 0];
C    = cell(1,3);
C{1} = spm_bsplinc(res.Twarp(:,:,:,1),prm);
C{2} = spm_bsplinc(res.Twarp(:,:,:,2),prm);
C{3} = spm_bsplinc(res.Twarp(:,:,:,3),prm);

mg  = res.mg;
mn  = res.mn;
vr  = res.vr;

do_defs = do_cls | do_defs;
if true,
    [pth,nam,ext1]=fileparts(res.image(1).fname);
    Ndef      = nifti;
    Ndef.dat  = file_array(fullfile(pth,['iy_', nam1, '.nii']),...
                                      [res.image(1).dim(1:3),1,3],...
                                      [spm_type('float32') spm_platform('bigend')],...
                                      0,1,0);
    Ndef.mat  = res.image(1).mat;
    Ndef.mat0 = res.image(1).mat;
    Ndef.descrip = ['Deformation'];
    create(Ndef);
end

spm_progress_bar('init',length(x3),['Working on ' nam],'Planes completed');
M = tpm.M\res.Affine*res.image(1).mat;

for z=1:length(x3),

    % Bias corrected image
    cr = cell(1,N);
    for n=1:N,
        f          = spm_sample_vol(res.image(n),x1,x2,o*x3(z),0);
        bf         = exp(transf(chan(n).B1,chan(n).B2,chan(n).B3(z,:),chan(n).T));
        cr{n}      = bf.*f;
        if ~isempty(chan(n).Nc),
            % Write a plane of bias corrected data
            chan(n).Nc.dat(:,:,z,chan(n).ind(1),chan(n).ind(2)) = cr{n};
        end;
        if ~isempty(chan(n).Nf),
            % Write a plane of bias field
            chan(n).Nf.dat(:,:,z,chan(n).ind(1),chan(n).ind(2)) = bf;
        end;
    end


    if do_defs,
        [t1,t2,t3] = defs(C,z,res.MT,prm,x1,x2,x3,M);
        if exist('Ndef','var'),
            M1 = tpm.M;
            tmp = M1(1,1)*t1 + M1(1,2)*t2 + M1(1,3)*t3 + M1(1,4);
            Ndef.dat(:,:,z,1,1) = tmp;
            tmp = M1(2,1)*t1 + M1(2,2)*t2 + M1(2,3)*t3 + M1(2,4);
            Ndef.dat(:,:,z,1,2) = tmp;
            tmp = M1(3,1)*t1 + M1(3,2)*t2 + M1(3,3)*t3 + M1(3,4);
            Ndef.dat(:,:,z,1,3) = tmp;
        end

        if exist('y','var'),
            y(:,:,z,1) = t1;
            y(:,:,z,2) = t2;
            y(:,:,z,3) = t3;
        end

        if do_cls,
            msk        = (f==0) | ~isfinite(f);
            q1  = likelihoods(cr,[],mg,mn,vr);
            q1  = reshape(q1,[d(1:2),numel(mg)]);

            q  = zeros([d(1:2) Kb]);
            b  = spm_sample_priors8(tpm,t1,t2,t3);
            for k1=1:Kb,
                q(:,:,k1) = sum(q1(:,:,find(lkp==k1)),3).*b{k1};
            end

            sq = sum(q,3);
            for k1=1:Kb,
                tmp            = q(:,:,k1);
                tmp(msk)       = 0;
                tmp            = tmp./sq;
                if ~isempty(cls{k1}),
                    cls{k1}(:,:,z) = uint8(round(255 * tmp));
                end
                if ~isempty(tiss(k1).Nt),
                    tiss(k1).Nt.dat(:,:,z,ind(1),ind(2)) = tmp;
                end
            end
        end
    end
    spm_progress_bar('set',z);
end
spm_progress_bar('clear');

NY = nifti(PY);
NT = nifti(PT);

y0 = NY.dat;
d0 = size(y0);d0 = d0(1:3);
M0 = NY.mat;

d1 = size(NT.dat); d1 = d1(1:3);
M1 = NT.mat;

% y maps from individual to mm coordinates in the template.  Change
% this so that it maps to voxel indices in the template instead.
y      = affind(squeeze(single(y0(:,:,:,:,:))),inv(M1));

% Generate a grid of mm coordinates for each voxel within the
% individual's image.
x      = affind(rgrid(d0(1:3)),M0);

% These mm coordinates are pushed to their new locations according
% to the original mapping (y).  Note that the resulting y is scaled
% at each point by the number of voxels that are mapped (w).
[y,w]  = dartel3('push',x,y,d1(1:3));

% Generate another grid of mm indices at each voxel in the template.
x      = affind(rgrid(d1(1:3)),M1);

% Fit a the closest affine transform through the mapping.
M      = getaff(x,y,w);

% Multiply the mm indices by this affine transform.  This transform
% is subtracted from the nonlinear warp so that a smooth displacement
% field can be determined.  The transform is then re-added to the
% displacement field.
x      = affind(x,M);
vx     = sqrt(sum(M1(1:3,1:3).^2));
for m=1:3,
    % Essentially, divide (y-x*w) by w but avoid divisions by zero by
    % biasing the result to be spatially smooth.
    y(:,:,:,m) = optimNn(w,y(:,:,:,m)-x(:,:,:,m).*w,[1  vx  0.001 0 0  2 1]) + x(:,:,:,m);
end




N      = nifti;
N.dat  = file_array('iy_junk.nii',[d1,1,3],'float32-be',0,1,0);
N.mat  = M1;
N.mat0 = M1;
N.descrip = 'Inverse Warp';
create(N);
N.dat(:,:,:,:,:) = reshape(y,[d1,1,3]);


P  = 'c1sMT04-0003-00001-000176-01_A.nii';
NC = nifti(P);
c  = single(NC.dat(:,:,:));

y      = affind(squeeze(single(y0(:,:,:,:,:))),inv(M1));
c  = single(NC.dat(:,:,:));

[c,w]  = dartel3('push',c,y,d1(1:3));
c      = optimNn(w,c,[2  vx  100 0 0  2 2]);

imagesc(c(:,:,35)'); axis image xy


return;
%=======================================================================

%=======================================================================
function [x1,y1,z1] = defs(sol,z,MT,prm,x0,y0,z0,M)
iMT = inv(MT);
x1  = x0*iMT(1,1)+iMT(1,4);
y1  = y0*iMT(2,2)+iMT(2,4);
z1  = (z0(z)*iMT(3,3)+iMT(3,4))*ones(size(x1));
x1a = x0    + spm_bsplins(sol{1},x1,y1,z1,prm);
y1a = y0    + spm_bsplins(sol{2},x1,y1,z1,prm);
z1a = z0(z) + spm_bsplins(sol{3},x1,y1,z1,prm);
x1  = M(1,1)*x1a + M(1,2)*y1a + M(1,3)*z1a + M(1,4);
y1  = M(2,1)*x1a + M(2,2)*y1a + M(2,3)*z1a + M(2,4);
z1  = M(3,1)*x1a + M(3,2)*y1a + M(3,3)*z1a + M(3,4);
return;
%=======================================================================

%=======================================================================
function t = transf(B1,B2,B3,T)
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t = zeros(size(B1,1),size(B2,1),size(B3,1));
end;
return;
%=======================================================================

%=======================================================================
function dat = decimate(dat,fwhm)
% Convolve the volume in memory (fwhm in voxels).
lim = ceil(2*fwhm);
x  = -lim(1):lim(1); x = spm_smoothkern(fwhm(1),x); x  = x/sum(x);
y  = -lim(2):lim(2); y = spm_smoothkern(fwhm(2),y); y  = y/sum(y);
z  = -lim(3):lim(3); z = spm_smoothkern(fwhm(3),z); z  = z/sum(z);
i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
spm_conv_vol(dat,dat,x,y,z,-[i j k]);
return;
%=======================================================================

%=======================================================================
function p = likelihoods(f,bf,mg,mn,vr)
K  = numel(mg);
N  = numel(f);
M  = numel(f{1});
cr = zeros(M,N);
for n=1:N,
    if isempty(bf),
        cr(:,n) = double(f{n}(:));
    else
        cr(:,n) = double(f{n}(:).*bf{n}(:));
    end
end
p  = ones(numel(f{1}),K);
for k=1:K,
    amp    = mg(k)/sqrt((2*pi)^N * det(vr(:,:,k)));
    d      = cr - repmat(mn(:,k)',M,1);
    p(:,k) = amp * exp(-0.5* sum(d.*(d/vr(:,:,k)),2));
end
p = p + 1024*eps;
%=======================================================================

%=======================================================================
function y1 = affind(y0,M);
y1 = zeros(size(y0),'single');
for d=1:3,
    y1(:,:,:,d) = y0(:,:,:,1)*M(d,1) + y0(:,:,:,2)*M(d,2) + y0(:,:,:,3)*M(d,3) + M(d,4);
end
%=======================================================================

%=======================================================================
function y1 = reweight(y0,w)
y1 = zeros(size(y0),'single');
for d=1:size(y0,4),
    y1(:,:,:,d) = y0(:,:,:,d).*w;
end
%=======================================================================

%=======================================================================
function x = rgrid(d);
x       = zeros([d(1:3),3],'single');
[x1,x2] = ndgrid(1:d(1),1:d(2));
for k=1:d(3),
    x(:,:,k,1) = x1;
    x(:,:,k,2) = x2;
    x(:,:,k,3) = k;
end
%=======================================================================

%=======================================================================
function [M,R] = getaff(x,y,w1,w2)
% Determine the affine (and rigid) transform mapping x to y
% FORMAT [M,R] = getaff(X,Y,W1,W2)
% X  - n1*n2*n3*3 array of floats representing coordinates.
% Y  - n1*n2*n3*3 array of floats representing coordinates.
% W1 - n1*n2*n3   array of floats representing weights.
% W2 - n1*n2*n3   array of floats representing weights.
%
% M  - an affine transform
% R  - a rigid-body transform
%
% The code treats X and Y as reshaped versions (n1*n2*n3) x 3,
% and W1 and W2 as column vectors.
%
% It generates XX = [diag(W1)*X W1]'*diag(W2)*[diag(W1)*X W1]
% and          XY = [diag(W1)*X W1]'*diag(W2)*[Y W1]
%
% These can then be used to compute an affine transform (M),
% by M = (XX\XY)'
% A weighted procrustes decomposition is also performed,
% so that a rigid-body transform matrix (R) is returned.
%
% If W1 or W2 are empty or not passed, then they are assumed
% to be all ones.

XX = zeros(4);
XY = zeros(4);
d  = size(x);
o  = ones(d(1)*d(2),1);
for k=1:size(x,3),
    xk  = reshape(x(:,:,k,:),[d(1)*d(2),3]);
    if nargin<3 && isempty(w1) && isempty(w2),
        ox = o;
        oy = o;
    else
        if nargin>=4 && ~isempty(w1) && ~isempty(w2),
            oy = reshape(wx(:,:,k), [d(1)*d(2),1]);
            ox = reshape( w(:,:,k), [d(1)*d(2),1]).*oy;
        elseif nargin>=3 && ~isempty(w1),
            ox = reshape(wx(:,:,k), [d(1)*d(2),1]);
            oy = ox;
        elseif nargin>=4 && ~isempty(w2),
            ox = reshape(wx(:,:,k), [d(1)*d(2),1]);
        end
        xk(:,1) = xk(:,1).*ox;
        xk(:,2) = xk(:,2).*ox;
        xk(:,3) = xk(:,3).*ox;
    end
    yk  = reshape(y(:,:,k,:),[d(1)*d(2),3]);
    msk = find(all(isfinite(xk),2) & all(isfinite(yk),2));
    X   = [xk(msk,:), ox(msk)];
    Y   = [yk(msk,:), oy(msk)];
    XX  = XX + double(X'*X);
    XY  = XY + double(X'*Y);
end
M = (XX\XY)';

% Procrustes decomposition
XX1 = XX - XX(:,4)*XX(:,4)'/XX(4,4);
XY1 = XY - XY(:,4)*XY(4,:) /XY(4,4);
Z   = (XX1(1:3,1:3)\XY1(1:3,1:3))';
[U,S,V] = svd(Z); % Decompose into rotate, zoom and rotate.
R   = [U*V' zeros(3,1);0 0 0 1];    % Rotation
T1  = [eye(4,3) -XY(:,4) /XY(4,4)]; % Initial translation of centre of mass to origin.
T2  = [eye(4,3) -XY(4,:)'/XY(4,4)]; % Final translation of origin to centre of mass.
R   = T2 * R * T1;

