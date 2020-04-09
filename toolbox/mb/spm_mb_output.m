function res = spm_mb_output(dat,mu,sett)
%__________________________________________________________________________
%
% Write output from groupwise normalisation and segmentation of images.
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% $Id$

if nargin==1
    load(dat);
    if isfield(sett.mu,'exist')
        mu = sett.mu.exist.mu;
    elseif isfield(sett.mu,'create')
        mu = sett.mu.create.mu;
    end
    mu = nifti(mu);
    mu = single(mu.dat(:,:,:,:,:));
end

% struct for saving paths of data written to disk
N   = numel(dat);
cl  = cell(N,1);
res = struct('inu',cl,'im',cl,'imc',cl,'c',cl,'wim',cl, ...
             'wimc',cl,'wc',cl,'mwc',cl);

for n=1:N % Loop over subjects
    res(n) = ProcessSubject(dat(n),res(n),mu,sett);
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
function resn = ProcessSubject(datn,resn,mu,sett)

% Parse function settings
dmu        = sett.mu.d;
Mmu        = sett.mu.Mmu;
dir_res    = sett.odir;
do_infer   = true;
mrf        = 1;
nit_mrf    = 8;
write_inu  = true; % field
write_im   = [true true true true]; % image, corrected, warped, warped corrected
write_tc   = [true true true]; % native, warped, warped-mod

% Get parameters
df     = datn.dm;
onam   = datn.onam;
Mn     = datn.Mat;
do_inu  = true;

if any(write_inu(:)) || any(write_im(:)) || any(write_tc(:)) 
    psi = spm_mb_io('GetData',datn.psi);
    psi = MatDefMul(psi,inv(Mmu));
end

if isfield(datn.model,'gmm') &&...
   (any(write_inu(:)) ||...
    any(write_im(:) ) ||...
    any(write_tc(:) ))
    % Input data were intensity images
    %------------------

    mun = spm_mb_shape('Pull1',mu,psi);
    mun = spm_mb_shape('TemplateK1',mun,4);
    mun = reshape(mun,size(mun,1)*size(mun,2)*size(mun,3),size(mun,4));

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

    % Get bias field modulated image data
    fn = inu.*fn;
    if do_infer
        % Infer missing values
        sample_post = do_infer > 1;
        MU          = gmm.m;
        A           = bsxfun(@times, gmm.W, reshape(gmm.n, [1 1 Kmg]));
        fn          = spm_gmm_lib('InferMissing',reshape(fn,[prod(df) C]),...
                                  zn,{MU,A},code_image,sample_post);
        clear code
    end

    % If using multiple Gaussians per tissue, collapse so that zn is of
    % size K1
    if Kmg > K1
        for k=1:K1,
            zn(:,k) = sum(zn(:,mg_ix==k),2);
        end
        zn(:,K1 + 1:end)    = [];
    end
    clear mun msk_chn

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
            spm_mb_io('WriteNii',fpth,inu(:,:,:,c),...
                      Mn,[descrip 'c=' num2str(c) ')']);
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
            spm_mb_io('WriteNii',fpth,fn(:,:,:,c)./inu(:,:,:,c),...
                      Mn,[descrip 'c=' num2str(c) ')'],'int16');
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
            spm_mb_io('WriteNii',fpth,fn(:,:,:,c),...
                      Mn,[descrip 'c=' num2str(c) ')'],'int16');
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
            spm_mb_io('WriteNii',fpth,zn(:,:,:,k),...
                      Mn,[descrip 'k=' num2str(k) ')'],'uint8');
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


if any(reshape(write_tc(:,[2 3]),[],1)) ||...
   any(reshape(write_im(:,[3 4]),[],1))
    % Write normalised images
    %------------------

    % For improved push - subsampling density in each dimension
    sd = spm_mb_shape('SampDens',Mmu,Mn);

    if isfield(datn.model,'gmm')
        disp('<- 1')
        if any(write_im(:,3))
            disp(['<-2 ' mum2str(C)])
            % Write normalised image
            descrip = 'Normalised image (';
            pths    = {};
            for c=1:C
                disp(c)
                if ~write_im(c,3), continue; end
                nam     = ['wim' num2str(c) '_' onam '.nii'];
                fpth    = fullfile(dir_res,nam);
                [img,cnt] = spm_mb_shape('Push1',fn(:,:,:,c)./inu(:,:,:,c),...
                                         psi,dmu,sd);
                spm_mb_io('WriteNii',fpth,img./(cnt + eps('single')),...
                          Mmu,[descrip 'c=' num2str(c) ')'],'int16');
                pths{end + 1} = fpth;
            end
            resn.wim = pths;
        end

        if any(write_im(:,4))
            % Write normalised image corrected
            descrip = 'Normalised image corrected (';
            pths    = {};
            for c=1:C
                if ~write_im(c,4), continue; end
                nam       = ['wimc' num2str(c) '_' onam '.nii'];
                fpth      = fullfile(dir_res,nam);
                [img,cnt] = spm_mb_shape('Push1',fn(:,:,:,c),psi,dmu,sd);
                spm_mb_io('WriteNii',fpth,img./(cnt + eps('single')),...
                          Mmu,[descrip 'c=' num2str(c) ')'],'int16');
                pths{end + 1} = fpth;
            end
            resn.wimc = pths;
        end
    end

    if any(write_tc(:,2)) || any(write_tc(:,3))
        pths    = {};
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
                    spm_mb_io('WriteNii',fpth,...
                        img./(cnt + eps('single')),...
                        Mmu,...
                        ['Norm. tiss. (k=' num2str(k) ')'],...
                        'uint8');
                end
                if write_tc(k,3)
                    % Write normalised modulated segmentation
                    kmwc = kmwc + 1;
                    fpth = fullfile(dir_res,['mwc' num2str(k) '_' onam '.nii']);
                    resn.mwc{kmwc} = fpth;
                    img  = img*abs(det(Mn(1:3,1:3))/det(Mmu(1:3,1:3)));
                    spm_mb_io('WriteNii',fpth,img,...
                              Mmu,...
                              ['Norm. mod. tiss. (k=' num2str(k) ')'],...
                              'int16');
                end
            end
        end
    end
end

end
%==========================================================================
function phi = MatDefMul(phi,M)
d   = size(phi);
phi = reshape(bsxfun(@plus,reshape(phi,[prod(d(1:3)),3])*M(1:3,1:3)',M(1:3,4)'),d);
end

