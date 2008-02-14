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
% $Id: spm_preproc_write8.m 1151 2008-02-14 17:36:47Z john $

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

do_defs = do_cls;
if true,
    [pth,nam,ext1]=fileparts(res.image(1).fname);
    for n=1:3,
        VY(n) = struct('fname',   fullfile(pth,['y_', nam1, '.nii']),...
                       'dim',     res.image(1).dim(1:3),...
                       'dt',      [spm_type('float32') spm_platform('bigend')],...
                       'pinfo',   [1 0 0]',...
                       'mat',     res.image(1).mat,...
                       'n',       [1 n],...
                       'descrip', 'Deformation');
    end
    VY = spm_create_vol(VY);
    do_defs = true;
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
        if true,
            M1 = tpm.M;
            tmp = M1(1,1)*t1 + M1(1,2)*t2 + M1(1,3)*t3 + M1(1,4);
            spm_write_plane(VY(1),tmp,z);
            tmp = M1(2,1)*t1 + M1(2,2)*t2 + M1(2,3)*t3 + M1(2,4);
            spm_write_plane(VY(2),tmp,z);
            tmp = M1(3,1)*t1 + M1(3,2)*t2 + M1(3,3)*t3 + M1(3,4);
            spm_write_plane(VY(3),tmp,z);
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
function [g,w,c] = clean_gwc(g,w,c, level)
if nargin<4, level = 1; end;

b    = w;
b(1) = w(1);

% Build a 3x3x3 seperable smoothing kernel
%-----------------------------------------------------------------------
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

th1 = 0.15;
if level==2, th1 = 0.2; end;
% Erosions and conditional dilations
%-----------------------------------------------------------------------
niter = 32;
spm_progress_bar('Init',niter,'Extracting Brain','Iterations completed');
for j=1:niter,
        if j>2, th=th1; else th=0.6; end; % Dilate after two its of erosion.
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

