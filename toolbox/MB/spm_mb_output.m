function res = spm_mb_output(cfg)
% Write output from groupwise normalisation and segmentation of images
% FORMAT res = spm_mb_output(cfg)
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


res  = cfg.result;
if iscell(res), res = res{1}; end
if ischar(res), res  = load(res); end
sett = res.sett;
dat  = res.dat;

if isfield(sett.mu,'exist')
    mu = sett.mu.exist.mu;
elseif isfield(sett.mu,'create')
    mu = sett.mu.create.mu;
end
if mu(1) == '.', mu = fullfile(fileparts(cfg.result{1}), mu); end
mu = nifti(mu);
mu = single(mu.dat(:,:,:,:,:));

% If SPM has been compiled with OpenMP support then the number of threads
% are here set to speed up the algorithm
%--------------------------------------------------------------------------
nw   = spm_mb_shape('get_num_workers',sett,max(27,sett.K*5+17));
if sett.nworker > 1
    setenv('SPM_NUM_THREADS',sprintf('%d',0));
else
    setenv('SPM_NUM_THREADS',sprintf('%d',-1));
end

% struct for saving paths of data written to disk
N   = numel(dat);
cl  = cell(N,1);
res = struct('inu',cl,'i',cl,'mi',cl,'c',cl,'wi',cl, ...
             'wmi',cl,'wc',cl,'mwc',cl,'sm',cl);

write_tc = false(sett.K+1,4);
ind = cfg.c;   ind = ind(ind>=1 & ind<=sett.K+1); write_tc(ind,1) = true;
ind = cfg.wc;  ind = ind(ind>=1 & ind<=sett.K+1); write_tc(ind,2) = true;
ind = cfg.mwc; ind = ind(ind>=1 & ind<=sett.K+1); write_tc(ind,3) = true;
ind = cfg.sm;  ind = ind(ind>=1 & ind<=sett.K+1); write_tc(ind,4) = true;

opt = struct('write_inu',cfg.inu,...
             'write_im',[cfg.i cfg.mi cfg.wi cfg.wmi],...
             'write_tc',write_tc,...
             'mrf',cfg.mrf,...
             'vx',cfg.vox,...
             'bb',cfg.bb,...
             'odir',cfg.odir,...
             'fwhm',cfg.fwhm);
opt.proc_zn = cfg.proc_zn;

if nw > 1 && numel(dat) > 1 % PARFOR
    fprintf('Write output: ');
    parfor(n=1:N,nw)
        fprintf('.');
        res(n) = ProcessSubject(dat(n),mu,sett,opt);
    end
    fprintf(' done!\n');
else
    spm_progress_bar('Init',N,'Writing MB output','Subjects complete');
    for n=1:N % FOR
        res(n) = ProcessSubject(dat(n),mu,sett,opt);
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
if nargin < 4, nit = 10; end
P   = zeros(size(zn),'uint8');
G   = ones([size(zn,4),1],'single')*strength;
vx  = sqrt(sum(Mn(1:3,1:3).^2));
vx2 = 1./single(vx);
for i=1:nit
    spm_mrf(P,zn,G,vx2);
end
zn = single(P)/255;
%==========================================================================

%==========================================================================
% ProcessSubject()
function resn = ProcessSubject(datn,mu,sett,opt)

% Parse function settings
dmu         = sett.ms.d;
Mmu         = sett.ms.Mmu;
if isempty(opt.odir)
    dir_res = sett.odir;
else
    dir_res = opt.odir;
    if ~(exist(dir_res,'dir') == 7)
        mkdir(dir_res)
    end
end
do_infer    = true;
mrf         = opt.mrf;
write_inu   = opt.write_inu; % field
write_im    = opt.write_im;  % image, corrected, warped, warped corrected
write_tc    = opt.write_tc;  % native, warped, warped-mod, scalar momentum
fwhm        = opt.fwhm;      % FWHM for smoothing of warped tissues
vx          = opt.vx;        % Template space voxel size
bb          = opt.bb;        % Template space bounding box
proc_zn     = opt.proc_zn;   % Function for processing native space responsibilities

cl   = cell(1,1);
resn = struct('inu',cl,'i',cl,'mi',cl,'c',cl,'wi',cl, ...
              'wmi',cl,'wc',cl,'mwc',cl,'sm',cl);

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
        chan = spm_mb_appearance('inu_basis',gmm.T,df,datn.Mat);
        inu  = spm_mb_appearance('inu_field',gmm.T,chan);
        clear chan

        if any(write_inu == true) && any(do_inu == true)
            % Write bias field
            inu      = reshape(inu,[df(1:3) C]);
            resn.inu = cell(1,sum(write_inu));
            c1       = 0;
            for c=1:C
                if ~write_inu(c,1), continue; end
                c1   = c1 + 1;
                nam  = sprintf('inu%d_%s.nii',c,onam);
                fpth = fullfile(dir_res,nam);
                write_nii(fpth,inu(:,:,:,c), Mn, sprintf('INU (%d)',c));
                c1           = c1 + 1;
                resn.inu{c1} = fpth;
            end
        end
    else
        inu   = ones([1 C],'single');
    end
end

if any(write_im(:)) || any(write_tc(:))
    Mat = Mmu\spm_dexpm(double(datn.q),sett.B)*datn.Mat;
    psi = spm_mb_io('get_data',datn.psi);
    psi = MatDefMul(psi,inv(Mmu));
    psi = spm_mb_shape('compose',psi,spm_mb_shape('affine',datn.dm,Mat));
end


if isfield(datn.model,'gmm') && (any(write_im(:)) || any(write_tc(:)))

    % Get image(s)
    fn     = spm_mb_io('get_image',gmm);

    % Get warped tissue priors
    mun    = spm_mb_shape('pull1',mu,psi);
    mun    = spm_mb_classes('template_k1',mun,datn.delta);

    % Format for spm_gmm
    chan                     = spm_mb_appearance('inu_basis',gmm.T,df,datn.Mat,ones(1,C));
    [~,mf,vf]                = spm_mb_appearance('inu_recon',fn,[],chan,gmm.T,gmm.Sig);

    % Get label data
    if isa(datn.lab,'struct')
        label = uint8(spm_mb_io('get_data', datn.lab.f));
    else
        label = [];
    end
    lnP = dirichlet_logexpect(gmm.Alpha);

    % For improved push - subsampling density in each dimension
    sd  = spm_mb_shape('samp_dens',Mmu,Mn);

    zn  = spm_gmmlib('resp',gmm.m,gmm.b,gmm.W,gmm.nu,gmm.gam,...
                     uint64(mg_ix), mun,mf,vf, uint64(gmm.samp), label,lnP);
    zn  = cat(4,zn,1-sum(zn,4));

    if mrf > 0
        % Ad-hoc MRF clean-up of segmentation
        zn = PostProcMRF(zn,Mn,mrf);
    end

    if iscell(proc_zn) && ~isempty(proc_zn) && isa(proc_zn{1},'function_handle')
        % Applies a function that processes the native space responsibilities
        try
            zn = proc_zn{1}(zn);
        catch
            warning('Incorrect definition of out.proc_zn, no processing performed.')
        end
    end

    if any(write_tc(:,1) == true)
        % Write segmentations
        resn.c  = cell(1,sum(write_tc(:,1)));
        k1      = 0;
        for k=1:K1
            if ~write_tc(k,1), continue; end
            nam  = sprintf('c%.2d_%s.nii',k,onam);
            fpth = fullfile(dir_res,nam);
            write_nii(fpth,zn(:,:,:,k), Mn, sprintf('Tissue (%d)',k), 'uint8');
            k1         = k1 + 1;
            resn.c{k1} = fpth;
        end
    end


    if do_infer
        % Infer missing values
        mf = spm_gmmlib('infer', gmm.m,gmm.b,gmm.W,gmm.nu,gmm.gam, uint64(mg_ix), ...
                                 mun,mf,vf, uint64([1 1 1]), label,lnP);
    end

    if any(write_im(:,1))
        % Write image
        resn.i = cell(1,sum(write_im(:,1)));
        c1     = 0;
        for c=1:C
            if ~write_im(c,1), continue; end
            nam  = sprintf('i%d_%s.nii',c,onam);
            fpth = fullfile(dir_res,nam);
            write_nii(fpth,mf(:,:,:,c)./inu(:,:,:,c), Mn, sprintf('Image (%d)',c), 'int16');
            c1         = c1 + 1;
            resn.i{c1} = fpth;
        end
    end

    if any(write_im(:,2))
        % Write image corrected
        resn.mi = cell(1,sum(write_im(:,2)));
        c1      = 0;
        for c=1:C
            if ~write_im(c,2), continue; end
            nam  = sprintf('mi%d_%s.nii',c,onam);
            fpth = fullfile(dir_res,nam);
            write_nii(fpth, mf(:,:,:,c), Mn, sprintf('INU corr. (%d)',c), 'int16');
            c1          = c1 + 1;
            resn.mi{c1} = fpth;
        end
    end

    if any(write_im(:,3))
        % Write normalised image
        resn.wi = cell(1,sum(write_im(:,3)));
        c1      = 0;
        for c=1:C
            if ~write_im(c,3), continue; end
            nam       = sprintf('wi%d_%s.nii',c,onam);
            fpth      = fullfile(dir_res,nam);
            [img,cnt] = spm_mb_shape('push1',mf(:,:,:,c)./inu(:,:,:,c), psi,dmu,sd);
            write_nii(fpth,img./(cnt + eps('single')), Mmu, sprintf('Norm. (%d)',c), 'int16');
            clear img cnt
            c1           = c1 + 1;
            resn.wi{c1} = fpth;
        end
    end
    clear inu

    if any(write_im(:,4))
        % Write normalised image corrected
        resn.wmi = cell(1,sum(write_im(:,4)));
        c1       = 0;
        for c=1:C
            if ~write_im(c,4), continue; end
            nam       = sprintf('wmi%d_%s.nii',c,onam);
            fpth      = fullfile(dir_res,nam);
            [img,cnt] = spm_mb_shape('push1',mf(:,:,:,c),psi,dmu,sd);
            write_nii(fpth,img./(cnt + eps('single')), Mmu, sprintf('Norm. INU corr. (%d)',c),'int16');
            clear img cnt
            c1           = c1 + 1;
            resn.wmi{c1} = fpth;
        end
    end
    clear mf vf
end

if isfield(datn.model,'cat') && (any(write_tc(:,2)) || any(write_tc(:,3)))
    % Input data were segmentations
    %------------------
    zn = spm_mb_io('get_data',datn.model.cat.f);
    zn = cat(4,zn,1 - sum(zn,4));
    K1 = sett.K+1;
end

% For improved push - subsampling density in each dimension
sd    = spm_mb_shape('samp_dens',Mmu,Mn);
vx_mu = sqrt(sum(Mmu(1:3,1:3).^2));

if any(write_tc(:,2)) || any(write_tc(:,3)) || any(write_tc(:,4))
    if any(write_tc(:,2)), resn.wc  = cell(1,sum(write_tc(:,2))); end
    if any(write_tc(:,3)), resn.mwc = cell(1,sum(write_tc(:,3))); end
    if any(write_tc(:,4)), resn.sm  = cell(1,sum(write_tc(:,4))); end
    kwc  = 0;
    kmwc = 0;
    ksm  = 0;
    if any(write_tc(:,4))
        % The scalar momentum residuals are computed in native space and then
        % pushed. We therefore here compute the softmaxed K + 1 template in
        % native space.
        mun = spm_mb_shape('pull1',mu,psi);
        clear mu
        mun = spm_mb_shape('softmax0',mun,4);
        mun = cat(4,mun,max(1 - sum(mun,4),0));
    end

    % Possibly modify template space images' voxel size and FOV
    [Mmu,dmu,vx_mu,psi] = modify_fov(bb,vx,Mmu,dmu,vx_mu,psi,sett);

    % Write output
    for k=1:K1
        if write_tc(k,2) || write_tc(k,3) || write_tc(k,4)
            if write_tc(k,2) || write_tc(k,3)
                [img,cnt] = spm_mb_shape('push1',zn(:,:,:,k),psi,dmu,sd);
            end
            if write_tc(k,2)
                % Write normalised segmentation
                kwc          = kwc + 1;
                fpth         = fullfile(dir_res, sprintf('wc%.2d_%s.nii',k,onam));
                resn.wc{kwc} = fpth;
                wimg         = img./(cnt + eps('single'));
                spm_smooth(wimg,wimg,fwhm(1)./vx_mu);  % Smooth
                write_nii(fpth, wimg, Mmu, sprintf('Norm. tissue (%d)',k), 'uint8');
            end
            if write_tc(k,3)
                % Write normalised modulated segmentation
                kmwc           = kmwc + 1;
                fpth           = fullfile(dir_res,sprintf('mwc%.2d_%s.nii',k,onam));
                resn.mwc{kmwc} = fpth;
                wimg           = img*abs(det(Mn(1:3,1:3))/det(Mmu(1:3,1:3)));
                spm_smooth(wimg,wimg,fwhm(2)./vx_mu);  % Smooth
                write_nii(fpth,wimg, Mmu, sprintf('Norm. mod. tissue (%d)',k), 'int16');
            end
            clear img cnt
            if write_tc(k,4)
                % Write scalar momentum, reference:
                % "A comparison of various MRI feature types for characterizing
                %  whole brain anatomical differences using linear pattern
                %  recognition methods." Monte-Rubio, et al. NeuroImage (2018)
                ksm          = ksm + 1;
                fpth         = fullfile(dir_res,sprintf('sm%.2d_%s.nii',k,onam));
                resn.sm{ksm} = fpth;
                % Compute scalar momentum
                wimg                   = spm_mb_shape('push1',zn(:,:,:,k) - mun(:,:,:,k),psi,dmu,sd);
                wimg(~isfinite(wimg))  = 0;           % Assume all values are zero outside FOV
                spm_smooth(wimg,wimg,fwhm(3)./vx_mu); % Smooth
                write_nii(fpth,wimg, Mmu, sprintf('Scalar momentum (%d)',k), 'int16');
            end
            clear wimg
        end
    end
end
%==========================================================================

%==========================================================================
function lnP = dirichlet_logexpect(Alpha)
% Expectation of parameter lobs from Dirichlet distribution
% Note: this is a separate function because there is a variable psi
lnP = bsxfun(@minus, psi(Alpha), psi(sum(Alpha,1)));
if ~isempty(Alpha), lnP(1,:) = 0; end
%==========================================================================

%==========================================================================
function [Mmu,dmu,vx_mu,psi] = modify_fov(bb_out,vx_out,Mmu,dmu,vx_mu,psi,sett)
if any(isfinite(bb_out(:))) || any(isfinite(vx_out))
    % Get bounding-box
    if isfield(sett.mu,'exist')
        [bb0,vox0] = spm_get_bbox(sett.mu.exist.mu,  'old');
    else
        [bb0,vox0] = spm_get_bbox(sett.mu.create.mu, 'old');
    end
    vx_out     = vx_out(1)*ones(1,3);
    msk    = ~isfinite(vx_out); vx_out(msk) = vox0(msk);
    msk    = ~isfinite(bb_out); bb_out(msk) =  bb0(msk);
    bb_out = sort(bb_out);
    vx_out = abs(vx_out);
    % Adjust bounding box slightly - so it rounds to closest voxel.
    bb_out(:,1) = round(bb_out(:,1)/vx_out(1))*vx_out(1);
    bb_out(:,2) = round(bb_out(:,2)/vx_out(2))*vx_out(2);
    bb_out(:,3) = round(bb_out(:,3)/vx_out(3))*vx_out(3);
    dim = round(diff(bb_out)./vx_out+1);
    of  = -vx_out.*(round(-bb_out(1,:)./vx_out)+1);
    mat = [vx_out(1) 0 0 of(1) ; 0 vx_out(2) 0 of(2) ; 0 0 vx_out(3) of(3) ; 0 0 0 1];
    if det(Mmu(1:3,1:3)) < 0
        mat = mat*[-1 0 0 dim(1)+1; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    end
    M0    = mat\Mmu;
    Mmu   = mat;
    vx_mu = sqrt(sum(Mmu(1:3,1:3).^2));
    dmu   = dim;
    % Modify deformation
    psi = MatDefMul(psi,M0);
end
%==========================================================================

%==========================================================================
function phi = MatDefMul(phi,M)
d   = size(phi);
phi = reshape(bsxfun(@plus,reshape(phi,[prod(d(1:3)),3])*M(1:3,1:3)',M(1:3,4)'),d);
%==========================================================================

%==========================================================================
function write_nii(f,img,M,descrip,typ)
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
    error(sprintf('Can''t do datatype "%s"', typ));
end
Nii         = nifti;
Nii.dat     = fa;
Nii.mat     = M;
Nii.mat0    = M;
Nii.descrip = descrip;
create(Nii);
Nii.dat(:,:,:,:,:,:) = img;
%==========================================================================
