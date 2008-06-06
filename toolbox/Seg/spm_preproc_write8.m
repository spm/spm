function cls = spm_preproc_write8(res,tc,bf,df)
% Write out VBM preprocessed data
% FORMAT cls = spm_preproc_write8(res,tc,bf,df)
%____________________________________________________________________________
% Copyright (C) 2008 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: spm_preproc_write8.m 1798 2008-06-06 16:25:56Z john $

tpm = res.tpm;
if ~isstruct(tpm) || ~isfield(tpm, 'bg'),
    tpm = spm_load_priors8(tpm);
end

lkp = res.lkp;
Kb  = max(lkp);
N   = numel(res.image);
if nargin<2, tc = true(Kb,4); end % native, import, warped, warped-mod
if nargin<3, bf = true(N,2);  end % corrected, field
if nargin<4, df = true(1,2);  end % inverse, forward

[pth,nam]=fileparts(res.image(1).fname);
ind  = res.image(1).n;
d    = res.image(1).dim(1:3);

[x1,x2,o] = ndgrid(1:d(1),1:d(2),1);
x3  = 1:d(3);

chan(N) = struct('B1',[],'B2',[],'B3',[],'T',[],'Nc',[],'Nf',[],'ind',[]);
for n=1:N,
    d3         = [size(res.Tbias{n}) 1];
    chan(n).B3 = spm_dctmtx(d(3),d3(3),x3);
    chan(n).B2 = spm_dctmtx(d(2),d3(2),x2(1,:)');
    chan(n).B1 = spm_dctmtx(d(1),d3(1),x1(:,1));
    chan(n).T  = res.Tbias{n};

    [pth1,nam1,ext1] = fileparts(res.image(n).fname);
    chan(n).ind      = res.image(n).n;

    if bf(n,1),
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

    if bf(n,2),
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

do_cls   = any(tc(:)) || nargout>1;
tiss(Kb) = struct('Nt',[]);
cls      = cell(1,Kb);
for k1=1:Kb,
    if tc(k1,4) || any(tc(:,3)) || nargout>=1,
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

do_defs = any(df);
do_defs = do_cls | do_defs;
if do_defs,
    if df(2),
        [pth,nam,ext1]=fileparts(res.image(1).fname);
        Ndef      = nifti;
        Ndef.dat  = file_array(fullfile(pth,['iy_', nam1, '.nii']),...
                               [res.image(1).dim(1:3),1,3],...
                               [spm_type('float32') spm_platform('bigend')],...
                               0,1,0);
        Ndef.mat  = res.image(1).mat;
        Ndef.mat0 = res.image(1).mat;
        Ndef.descrip = 'Inverse Deformation';
        create(Ndef);
    end
    if df(1) || any(tc(:,3)) || any(tc(:,4)) || nargout>=1,
        y = zeros([res.image(1).dim(1:3),3],'single');
    end
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
            msk = (f==0) | ~isfinite(f);
            q1  = likelihoods(cr,[],mg,mn,vr);
            q1  = reshape(q1,[d(1:2),numel(mg)]);

            q   = zeros([d(1:2) Kb]);
            b   = spm_sample_priors8(tpm,t1,t2,t3);
            for k1=1:Kb,
                q(:,:,k1) = sum(q1(:,:,lkp==k1),3).*b{k1};
            end

            sq = sum(q,3) + eps^2;
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
    

M0 = res.image(1).mat;

d1 = size(tpm.dat{1}); d1 = d1(1:3);
M1 = tpm.M;

if any(tc(:,3)),
    C = zeros([d1,Kb],'single');
end

if any(tc(:,3)) || any(tc(:,4)) || nargout>=1,
    spm_progress_bar('init',Kb,'Warped Tissue Classes','Classes completed');
    for k1 = 1:Kb,
        if ~isempty(cls{k1}),
            c = single(cls{k1})/255;
            if any(tc(:,3)),
                [c,w]  = dartel3('push',c,y,d1(1:3));
            else
                c      = dartel3('push',c,y,d1(1:3));
            end
            if nargout>=1,
                cls{k1} = c;
            end
            if tc(k1,4),
                N      = nifti;
                N.dat  = file_array(fullfile(pth,['mwc', num2str(k1), nam, ext1]),...
                                    d1,...
                                    [spm_type('float32') spm_platform('bigend')],...
                                    0,1,0);
                N.mat  = M1;
                N.mat0 = M1;
                N.descrip = ['Jac. sc. warped tissue class ' num2str(k1)];
                create(N);
                N.dat(:,:,:) = c*abs(det(M0(1:3,1:3))/det(M1(1:3,1:3)));
            end
            if any(tc(:,3)),
                vx          = sqrt(sum(M1(1:3,1:3).^2));
                C(:,:,:,k1) = optimNn(w,c,[1  vx  1e-4 1e-6 0  3 2]);
            end
            spm_progress_bar('set',k1);
        end
    end
    spm_progress_bar('Clear');
end

if any(tc(:,3)),
    spm_progress_bar('init',Kb,'Writing Warped Tis Cls','Classes completed');
    C = max(C,eps);
    s = sum(C,4);
    for k1=1:Kb,
        if tc(k1,3),
            N      = nifti;
            N.dat  = file_array(fullfile(pth,['wc', num2str(k1), nam, ext1]),...
                                d1,'uint8-be',0,1/255,0);
            N.mat  = M1;
            N.mat0 = M1;
            N.descrip = ['Warped tissue class ' num2str(k1)];
            create(N);
            N.dat(:,:,:) = C(:,:,:,k1)./s;
        end
        spm_progress_bar('set',k1);
    end
    spm_progress_bar('Clear');
    clear C s
end

if df(1),
    y         = spm_invert_def(y,M1,d1,M0,[1 0]);
    N         = nifti;
    N.dat     = file_array(fullfile(pth,['y_', nam1, '.nii']),...
                           [d1,1,3],'float32-be',0,1,0);
    N.mat     = M1;
    N.mat0    = M1;
    N.descrip = 'Deformation';
    create(N);
    N.dat(:,:,:,:,:) = reshape(y,[d1,1,3]);
end

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
%=======================================================================

