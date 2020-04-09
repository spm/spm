function res = spm_mb_output(dat,mu,sett)
%__________________________________________________________________________
%
% Write output from groupwise normalisation and segmentation of images.
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% struct for saving paths of data written to disk
N   = numel(dat);
cl  = cell(N,1);
res = struct('inu',cl,'im',cl,'imc',cl,'c',cl,'y',cl,'wim',cl, ...
             'wimc',cl,'wc',cl,'mwc',cl,'v',cl);

for n=1:N % Loop over subjects
    res(n) = ProcessSubject(dat(n),res(n),mu,n,sett);
end
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
function resn = ProcessSubject(datn,resn,mu,ix,sett)

% Parse function settings
clean_def  = false;
dmu        = sett.ms.d;
Mmu        = sett.ms.Mmu;
dir_res    = sett.dir;
do_infer   = true;
mrf        = 0;
nit_mrf    = 8;
write_inu   = true; % field
write_df   = [true false]; % forward, inverse
write_im   = [true true true true]; % image, corrected, warped, warped corrected
write_tc   = [true true true]; % native, warped, warped-mod
write_vel  = false;

% Get parameters
df     = datn.dm;
onam   = datn.onam;
Mr     = spm_dexpm(double(datn.q),sett.B);
Mn     = datn.Mat;
do_inu  = true;

if ~(all(write_inu(:) == false) && all(write_im(:) == false) && all(write_tc(:) == false) && all(write_df(:) == false))
    psi0 = spm_mb_io('GetData',datn.psi);
end

psi = spm_mb_io('GetData',datn.psi);
psi = spm_mb_shape('Compose',psi,spm_mb_shape('Affine',df,inv(Mmu));


if isfield(datn.model,'gmm') && (any(write_inu(:) == true) || any(write_im(:) == true) || any(write_tc(:) == true))
    % Input data were intensity images
    %------------------

    mun = spm_mb_shape('Pull1',mu,psi);
    mun = spm_mb_shape('TemplateK1',mun,4);

    gmm    = datn.model.gmm;
    gmms   = sett.gmm(gmm.pop);
    C      = gmms.C;
    mg_ix  = gmms.mg_ix;
    K      = sett.K;
    K1     = K + 1;
    Kmg    = numel(mg_ix);

    % Integrate K1 and C into write settings
    if size(write_inu,1) == 1 && C  > 1, write_inu = repmat(write_inu,[C  1]); end
    if size(write_im, 1) == 1 && C  > 1, write_im  = repmat(write_im, [C  1]); end
    if size(write_tc, 1) == 1 && K1 > 1, write_tc  = repmat(write_tc, [K1 1]); end


fprintf('Computing.\n');

    if any(do_inu == true)
        % Get bias field
        chan = spm_mb_appearance('BiasBasis',gmm.T,df,datn.Mat,ones(1,C));
        inu  = spm_mb_appearance('BiasField',gmm.T,chan);
    else
        inu   = ones([1 C],'single');
    end

    % Get image(s)
    fn      = spm_mb_io('GetImage',gmm);

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

    % GMM posterior
    m = gmm.m;
    b = gmm.b;
    W = gmm.W;
    n = gmm.n;

    % Get responsibilities
    zn  = spm_mb_appearance('Responsibility',m,b,W,n,inufn,mun,msk_chn);
    zn  = spm_gmm_lib('cell2obs', zn, code_image, msk_chn);
    clear mun msk_chn

    % Get bias field modulated image data
    fn = inu.*fn;
    if do_infer
        % Infer missing values
        sample_post = do_infer > 1;
        MU          = gmm.m;
        A           = bsxfun(@times, gmm.W, reshape(gmm.n, [1 1 Kmg]));
        fn          = spm_gmm_lib('InferMissing',reshape(fn,[prod(df) C]),zn,{MU,A},code_image,sample_post);
        clear code
    end

    % If using multiple Gaussians per tissue, collapse so that zn is of
    % size K1
    if Kmg > K1
        for k=1:K1, zn(:,k) = sum(zn(:,mg_ix==k),2); end
        zn(:,K1 + 1:end)    = [];
    end

    % Make 3D
    if any(do_inu == true)
        inu = reshape(inu,[df(1:3) C]);
    else
        inu = reshape(inu,[1 1 1 C]);
    end
    fn = reshape(fn,[df(1:3) C]);
    zn = reshape(zn,[df(1:3) K1]);

    if mrf > 0
        % Ad-hoc MRF clean-up of segmentation
        zn = PostProcMRF(zn,Mn,mrf,nit_mrf);
    end

    if any(write_inu == true) && any(do_inu == true)
        % Write bias field
        descrip = 'Bias field (';
        pths    = {};
        for c=1:C
            if ~write_inu(c,1), continue; end
            nam  = ['inu' num2str(c) '_' onam '.nii'];
            fpth = fullfile(dir_res,nam);
            spm_mb_io('WriteNii',fpth,inu(:,:,:,c),Mn,[descrip 'c=' num2str(c) ')']);
            pths{end + 1} = fpth;
        end
        resn.inu = pths;
    end

    if any(write_im(:,1) == true)
        % Write image
        descrip = 'Image (';
        pths    = {};
        for c=1:C
            if ~write_im(c,1), continue; end
            nam  = ['im' num2str(c) '_' onam '.nii'];
            fpth = fullfile(dir_res,nam);
            spm_mb_io('WriteNii',fpth,fn(:,:,:,c)./inu(:,:,:,c),Mn,[descrip 'c=' num2str(c) ')'],'int16');
            pths{end + 1} = fpth;
        end
        resn.im = pths;

        % Write image corrected
        descrip = 'Image corrected (';
        pths    = {};
        for c=1:C
            if ~write_im(c,2), continue; end
            nam  = ['imc' num2str(c) '_' onam '.nii'];
            fpth = fullfile(dir_res,nam);
            spm_mb_io('WriteNii',fpth,fn(:,:,:,c),Mn,[descrip 'c=' num2str(c) ')'],'int16');
            pths{end + 1} = fpth;
        end
        resn.imc = pths;
    end

    if any(write_tc(:,1) == true)
        % Write segmentations
        descrip = 'Tissue (';
        pths    = {};
        for k=1:K1
            if ~write_tc(k,1), continue; end
            nam  = ['c' num2str(k) '_' onam '.nii'];
            fpth = fullfile(dir_res,nam);
            spm_mb_io('WriteNii',fpth,zn(:,:,:,k),Mn,[descrip 'k=' num2str(k) ')'],'uint8');
            pths{end + 1} = fpth;
        end
        resn.c = pths;
    end
elseif isfield(datn.model,'cat')
    % Input data were segmentations
    %------------------
    zn = spm_mb_io('GetData',datn.model.cat.f);
    zn = cat(4,zn,1 - sum(zn,4));
end

fprintf('More writing.\n');

if any(write_df == true) || any(reshape(write_tc(:,[2 3]),[],1) == true) ||  any(reshape(write_im(:,[3 4]),[],1) == true)
    % Write forward deformation and/or normalised images
    %------------------

    % For imporved push - subsampling density in each dimension
    sd = spm_mb_shape('SampDens',Mmu,Mn);

    % Get forward deformation (pulls template into subject space)
    psi = spm_mb_shape('Compose',psi0,spm_mb_shape('Affine',df,Mmu\Mr*Mn));
    clear psi0

    if isfield(datn,'gmm') && any(write_im(:,3) == true)
        % Write normalised image
        descrip = 'Normalised image (';
        pths    = {};
        for c=1:C
            if ~write_im(c,3), continue; end
            nam     = ['wim' num2str(c) '_' onam '.nii'];
            fpth    = fullfile(dir_res,nam);
            [img,cnt] = spm_mb_shape('Push1',fn(:,:,:,c)./inu(:,:,:,c),psi,dmu,sd);
            spm_mb_io('WriteNii',fpth,img./(cnt + eps('single')),Mmu,[descrip 'c=' num2str(c) ')'],'int16');
            pths{end + 1} = fpth;
        end
        resn.wim = pths;
    end

    if isfield(datn,'gmm') && any(write_im(:,4) == true)
        % Write normalised image corrected
        descrip = 'Normalised image corrected (';
        pths    = {};
        for c=1:C
            if ~write_im(c,4), continue; end
            nam       = ['wimc' num2str(c) '_' onam '.nii'];
            fpth      = fullfile(dir_res,nam);
            [img,cnt] = spm_mb_shape('Push1',fn(:,:,:,c),psi,dmu,sd);
            spm_mb_io('WriteNii',fpth,img./(cnt + eps('single')),Mmu,[descrip 'c=' num2str(c) ')'],'int16');
            pths{end + 1} = fpth;
        end
        resn.wimc = pths;
    end

    if any(write_tc(:,2) == true)
        % Write normalised segmentations
        descrip = 'Normalised tissue (';
        pths    = {};
        for k=1:K1
            if ~write_tc(k,2), continue; end
            nam       = ['wc' num2str(k) '_' onam '.nii'];
            fpth      = fullfile(dir_res,nam);
            [img,cnt] = spm_mb_shape('Push1',zn(:,:,:,k),psi,dmu,sd);
            spm_mb_io('WriteNii',fpth,img./(cnt + eps('single')),Mmu,[descrip 'k=' num2str(k) ')'],'uint8');
            pths{end + 1} = fpth;
        end
        resn.wc = pths;
    end

    if any(write_tc(:,3) == true)
        % Write normalised modulated segmentations
        descrip = 'Normalised modulated tissue (';
        pths    = {};
        for k=1:K1
            if ~write_tc(k,3), continue; end
            nam  = ['mwc' num2str(k) '_' onam '.nii'];
            fpth = fullfile(dir_res,nam);
            img  = spm_mb_shape('Push1',zn(:,:,:,k),psi,dmu);
            img  = img*abs(det(Mn(1:3,1:3))/det(Mmu(1:3,1:3)));
            spm_mb_io('WriteNii',fpth,img,Mmu,[descrip 'k=' num2str(k) ')'],'int16');
            pths{end + 1} = fpth;
        end
        resn.mwc = pths;
    end

    if write_df(1)
        % Write forward deformation  (pulls template into subject space, use push to go other way)
       %psi     = reshape(reshape(psi,[prod(df) 3])*Mmu(1:3,1:3)' + Mmu(1:3,4)',[df 1 3]);
        psi     = reshape(bsxfun(@plus, reshape(psi,[prod(df) 3])*Mmu(1:3,1:3)', Mmu(1:3,4)'),[df 1 3]);
        descrip = 'Forward deformation';
        nam     = ['y_' onam '.nii'];
        fpth    = fullfile(dir_res,nam);
        spm_mb_io('WriteNii',fpth,psi,Mn,descrip);
        resn.y  = fpth;
    end
end

% Clean-up
%if clean_def  && isa(datn.psi,'nifti') && isfile(datn.psi.dat.fname),          delete(datn.psi.dat.fname); end
end
%==========================================================================
