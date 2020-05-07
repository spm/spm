function res = spm_mb_output(cfg)
%__________________________________________________________________________
%
% Write output from groupwise normalisation and segmentation of images.
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% $Id$

load(char(cfg.result));
if isfield(sett.mu,'exist')
    mu = sett.mu.exist.mu;
elseif isfield(sett.mu,'create')
    mu = sett.mu.create.mu;
end
mu = nifti(mu);
mu = single(mu.dat(:,:,:,:,:));

% struct for saving paths of data written to disk
N   = numel(dat);
cl  = cell(N,1);
res = struct('inu',cl,'i',cl,'mi',cl,'c',cl,'wi',cl, ...
             'wmi',cl,'wc',cl,'mwc',cl);

write_tc = false(sett.K+1,3);
ind = cfg.c;   ind = ind(ind>=1 & ind<=sett.K+1); write_tc(ind,1) = true;
ind = cfg.wc;  ind = ind(ind>=1 & ind<=sett.K+1); write_tc(ind,2) = true;
ind = cfg.mwc; ind = ind(ind>=1 & ind<=sett.K+1); write_tc(ind,3) = true;

opt = struct('write_inu',cfg.inu,...
             'write_im',[cfg.i cfg.mi cfg.wi cfg.wmi],...
             'write_tc',write_tc);

spm_progress_bar('Init',N,'Writing MB output','Subjects complete');
for n=1:N % Loop over subjects
    res(n) = ProcessSubject(dat(n),res(n),mu,sett,opt);
    spm_progress_bar('Set',n);
end
spm_progress_bar('Clear');
end
%==========================================================================

%==========================================================================
%
% Utility functions
%
%==========================================================================

%==========================================================================
% PostProcMRF()
function zn = PostProcMRF(zn,Mn,strength,nit)
P   = zeros(size(zn),'uint8');
G   = ones([size(zn,4),1],'single')*strength;
vx  = sqrt(sum(Mn(1:3,1:3).^2));
vx2 = 1./single(vx);
for i=1:nit
    spm_mrf(P,zn,G,vx2);
end
zn = single(P)/255;
end
%==========================================================================

%==========================================================================
% ProcessSubject()
function resn = ProcessSubject(datn,resn,mu,sett,opt)

% Parse function settings
dmu        = sett.mu.d;
Mmu        = sett.mu.Mmu;
dir_res    = sett.odir;
do_infer   = true;
mrf        = 0;
nit_mrf    = 0;
write_inu  = opt.write_inu; % field
write_im   = opt.write_im;  % image, corrected, warped, warped corrected
write_tc   = opt.write_tc;  % native, warped, warped-mod

if ((~any(write_inu(:)) && ~any(write_im(:))) || ~isfield(datn.model,'gmm')) && ~any(write_tc(:))
    return;
end

% Get parameters
df     = datn.dm;
onam   = datn.onam;
Mn     = datn.Mat;
do_inu = true;

if isfield(datn.model,'gmm')

    % Input data were intensity images
    %------------------

    gmm    = datn.model.gmm;
    gmms   = sett.gmm(gmm.pop);
    C      = gmms.C;
    mg_ix  = gmms.mg_ix;
    K      = sett.K;
    K1     = K + 1;
    Kmg    = numel(mg_ix);

    % Integrate K1 and C into write settings
    if size(write_inu,1) == 1
        write_inu = repmat(write_inu,[C  1]);
    end
    if size(write_im, 1) == 1
        write_im  = repmat(write_im, [C  1]);
    end
    if size(write_tc, 1) == 1
        write_tc  = repmat(write_tc, [K1 1]);
    end

    if any(do_inu == true)
        % Get bias field
        chan = spm_mb_appearance('BiasBasis',gmm.T,df,datn.Mat,ones(1,C));
        inu  = spm_mb_appearance('BiasField',gmm.T,chan);
        clear chan

        if any(write_inu == true) && any(do_inu == true)
            % Write bias field
            descrip  = 'Bias field (';
            inu      = reshape(inu,[df(1:3) C]);
            resn.inu = cell(1,sum(write_inu));
            c1       = 0;
            for c=1:C
                if ~write_inu(c,1), continue; end
                c1   = c1 + 1;
                nam  = ['inu' num2str(c) '_' onam '.nii'];
                fpth = fullfile(dir_res,nam);
                WriteNii(fpth,inu(:,:,:,c), Mn,[descrip 'c=' num2str(c) ')']);
                c1           = c1 + 1;
                resn.inu{c1} = fpth;
            end
        end
    else
        inu   = ones([1 C],'single');
    end
end

if any(write_im(:)) || any(write_tc(:))
    psi    = spm_mb_io('GetData',datn.psi);
    psi    = MatDefMul(psi,inv(Mmu));
end


if isfield(datn.model,'gmm') && any(write_im(:)) || any(write_tc(:))

    % Get image(s)
    fn     = spm_mb_io('GetImage',gmm);

    % Get warped tissue priors
    mun    = spm_mb_shape('Pull1',mu,psi);
    mun    = spm_mb_shape('TemplateK1',mun,4);
    mun    = reshape(mun,size(mun,1)*size(mun,2)*size(mun,3),size(mun,4));

    % Get labels
    labels = spm_mb_appearance('GetLabels',datn,sett);
    mun    = mun + labels;
    clear labels

    % Integrate use of multiple Gaussians per tissue
    mg_w = gmm.mg_w;
    mun  = mun(:,mg_ix);
   %mun  = mun + log(mg_w);
    mun  = bsxfun(@plus, mun, log(mg_w));

    % Format for spm_gmm
    [inufn,code_image,msk_chn] = spm_gmm_lib('obs2cell', reshape(inu.*fn,[prod(df) C]));
    mun                        = spm_gmm_lib('obs2cell', mun, code_image, false);

    % Get responsibilities
    zn  = spm_mb_appearance('Responsibility',gmm.m,gmm.b,gmm.W,gmm.n,inufn,mun,msk_chn);
    zn  = spm_gmm_lib('cell2obs', zn, code_image, msk_chn);
    clear mun msk_chn inufn

    % Get bias field modulated image data
    fn = inu.*fn;
    if do_infer
        % Infer missing values
        sample_post = do_infer > 1;
        A           = bsxfun(@times, gmm.W, reshape(gmm.n, [1 1 Kmg]));
        fn          = spm_gmm_lib('InferMissing',reshape(fn,[prod(df) C]),...
                                  zn,{gmm.m,A},code_image,sample_post);
        clear code
    end
    clear code_image

    fn = reshape(fn,[df(1:3) C]);

    if any(write_im(:,1))
        % Write image
        descrip = 'Image (';
        resn.im = cell(1,sum(write_im(:,1)));
        c1      = 0;
        for c=1:C
            if ~write_im(c,1), continue; end
            nam  = ['i' num2str(c) '_' onam '.nii'];
            fpth = fullfile(dir_res,nam);
            WriteNii(fpth,fn(:,:,:,c)./inu(:,:,:,c), Mn, [descrip 'c=' num2str(c) ')'], 'int16');
            c1          = c1 + 1;
            resn.m{c1} = fpth;
        end
    end

    if any(write_im(:,2))
        % Write image corrected
        descrip = 'Image corrected (';
        resn.mi = cell(1,sum(write_im(:,2)));
        c1      = 0;
        for c=1:C
            if ~write_im(c,2), continue; end
            nam  = ['mi' num2str(c) '_' onam '.nii'];
            fpth = fullfile(dir_res,nam);
            WriteNii(fpth, fn(:,:,:,c), Mn, [descrip 'c=' num2str(c) ')'], 'int16');
            c1           = c1 + 1;
            resn.mi{c1} = fpth;
        end
    end

    % For improved push - subsampling density in each dimension
    sd = spm_mb_shape('SampDens',Mmu,Mn);

    if any(write_im(:,3))
        % Write normalised image
        descrip = 'Normalised image (';
        resn.wi = cell(1,sum(write_im(:,3)));
        c1      = 0;
        for c=1:C
            if ~write_im(c,3), continue; end
            nam       = ['wi' num2str(c) '_' onam '.nii'];
            fpth      = fullfile(dir_res,nam);
            [img,cnt] = spm_mb_shape('Push1',fn(:,:,:,c)./inu(:,:,:,c), psi,dmu,sd);
            WriteNii(fpth,img./(cnt + eps('single')), Mmu, [descrip 'c=' num2str(c) ')'], 'int16');
            clear img cnt
            c1           = c1 + 1;
            resn.wi{c1} = fpth;
        end
    end
    clear inu

    if any(write_im(:,4))
        % Write normalised image corrected
        descrip  = 'Normalised image corrected (';
        resn.wmi = cell(1,sum(write_im(:,4)));
        c1       = 0;
        for c=1:C
            if ~write_im(c,4), continue; end
            nam       = ['wmi' num2str(c) '_' onam '.nii'];
            fpth      = fullfile(dir_res,nam);
            [img,cnt] = spm_mb_shape('Push1',fn(:,:,:,c),psi,dmu,sd);
            WriteNii(fpth,img./(cnt + eps('single')), Mmu,[descrip 'c=' num2str(c) ')'],'int16');
            clear img cnt
            c1           = c1 + 1;
            resn.wmi{c1} = fpth;
        end
    end
    clear fn

    % If using multiple Gaussians per tissue, collapse so that zn is of size K1
    if Kmg > K1
        for k=1:K1,
            zn(:,k) = sum(zn(:,mg_ix==k),2);
        end
        zn(:,K1 + 1:end)    = [];
    end
    zn = reshape(zn,[df(1:3) K1]);

    if mrf > 0
        % Ad-hoc MRF clean-up of segmentation
        zn = PostProcMRF(zn,Mn,mrf,nit_mrf);
    end

    if any(write_tc(:,1) == true)
        % Write segmentations
        descrip = 'Tissue (';
        resn.c  = cell(1,sum(write_tc(:,1)));
        k1      = 0;
        for k=1:K1
            if ~write_tc(k,1), continue; end
            nam  = ['c' num2str(k) '_' onam '.nii'];
            fpth = fullfile(dir_res,nam);
            WriteNii(fpth,zn(:,:,:,k), Mn, [descrip 'k=' num2str(k) ')'], 'uint8');
            k1         = k1 + 1;
            resn.c{k1} = fpth;
        end
    end
end

if isfield(datn.model,'cat') && (any(write_tc(:,2)) || any(write_tc(:,3)))
    % Input data were segmentations
    %------------------
    zn = spm_mb_io('GetData',datn.model.cat.f);
    zn = cat(4,zn,1 - sum(zn,4));
end


% For improved push - subsampling density in each dimension
sd = spm_mb_shape('SampDens',Mmu,Mn);

if any(write_tc(:,2)) || any(write_tc(:,3))
    if any(write_tc(:,2)), resn.wc  = cell(1,sum(write_tc(:,2))); end
    if any(write_tc(:,3)), resn.mwc = cell(1,sum(write_tc(:,3))); end
    kwc  = 0;
    kmwc = 0;
    for k=1:K1
        if write_tc(k,2) || write_tc(k,3)
            [img,cnt] = spm_mb_shape('Push1',zn(:,:,:,k),psi,dmu,sd);
            if write_tc(k,2)
                % Write normalised segmentation
                kwc  = kwc + 1;
                fpth = fullfile(dir_res, ['wc' num2str(k) '_' onam '.nii']);
                resn.wc{kwc} = fpth;
                WriteNii(fpth, img./(cnt + eps('single')),...
                         Mmu, ['Norm. tiss. (k=' num2str(k) ')'], 'uint8');
            end
            if write_tc(k,3)
                % Write normalised modulated segmentation
                kmwc = kmwc + 1;
                fpth = fullfile(dir_res,['mwc' num2str(k) '_' onam '.nii']);
                resn.mwc{kmwc} = fpth;
                img  = img*abs(det(Mn(1:3,1:3))/det(Mmu(1:3,1:3)));
                WriteNii(fpth,img, Mmu, ['Norm. mod. tiss. (k=' num2str(k) ')'], 'int16');
            end
            clear img cnt
        end
    end
end

end
%==========================================================================

%==========================================================================
function phi = MatDefMul(phi,M)
d   = size(phi);
phi = reshape(bsxfun(@plus,reshape(phi,[prod(d(1:3)),3])*M(1:3,1:3)',M(1:3,4)'),d);
end
%==========================================================================

%==========================================================================
function WriteNii(f,img,M,descrip,typ)
if nargin<5, typ = 'float32'; end
switch typ
case 'float32'
    fa = file_array(f,size(img),typ,0);
case 'uint8'
    mx = max(img(isfinite(img(:))));
    fa = file_array(f,size(img),typ,0,mx/255,0);
case 'int16'
    mx = max(img(isfinite(img(:))));
    mn = min(img(isfinite(img(:))));
    s  = max(mx/32767,-mn/32768);
    fa = file_array(f,size(img),typ,0,s,0);
otherwise
    error('Can''t do datatype "%s"', typ);
end
Nii         = nifti;
Nii.dat     = fa;
Nii.mat     = M;
Nii.mat0    = M;
Nii.descrip = descrip;
create(Nii);
Nii.dat(:,:,:,:,:,:) = img;
end
%==========================================================================

