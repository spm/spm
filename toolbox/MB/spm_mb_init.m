function [dat,sett] = spm_mb_init(cfg)
% Initialisation of Multi-Brain data structures
% FORMAT [dat,sett] = spm_mb_init(cfg)
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


[dat,sett] = mb_init1(cfg);

% Done if there are no GMMs to fit
if sum(cellfun(@(c)isfield(c,'gmm'),{dat.model}))==0, return; end
[sett,dat] = random_init(sett,dat);
%==========================================================================

%==========================================================================
function [dat,sett] = mb_init1(cfg)
sett      = cfg;
sett.odir = sett.odir{1};
if (numel(sett.odir>=1) && (sett.odir(1) == '.')) || (numel(sett.odir>=2) && strcmp(sett.odir(1:2),['.' filesep]))
    sett.odir = fullfile(pwd,sett.odir);
end
if ~isempty(sett.odir) && ~(exist(sett.odir, 'dir') == 7)
    mkdir(sett.odir);
end
if isfield(sett.mu,'exist')
    fnam = sett.mu.exist{1};
    sett.mu.exist = struct('mu',fnam);
    f    = nifti(fnam);
    dmu  = size(f(1).dat,[1 2 3]);
    K    = size(f(1).dat,4);
    Mmu  = f(1).mat;
else
    dmu  = [0 0 0];
    Mmu  = eye(4);
    K    = cfg.mu.create.K;
    sett.mu.create    = rmfield(sett.mu.create,'K');
    sett.mu.create.mu = fullfile(sett.odir,['mu_' cfg.onam '.nii']);
end
sett.mu.d   = dmu;
sett.mu.Mmu = Mmu;
sett.K      = K;

% Affine Lie algebra basis functions
if ~isempty(sett.aff)
    B   = spm_mb_shape('affine_bases',sett.aff);
else
    B   = zeros([4 4 0]);
end
sett.B = B;
sett   = rmfield(sett,'aff');
sett   = rmfield(sett,'cat');
sett   = rmfield(sett,'gmm');


% Count the total number of subjects (N) and ensure that all
% channels of each population have the same number of subjects.
N = 0;
if numel(cfg.cat)>=1
    N = numel(cfg.cat{1});
end
for p=1:numel(cfg.gmm)
    Np = 0;
    if numel(cfg.gmm(p).chan)>=1
        Np = numel(cfg.gmm(p).chan(1).images);
        N  = N + Np;
    end
    for c=2:numel(cfg.gmm(p).chan)
        if numel(cfg.gmm(p).chan(c).images)~=Np
            error(sprintf('Incompatible numbers of scans over channels (pop-%d).',p));
        end
    end
    if isfield(cfg.gmm(p).labels,'true')
        if numel(cfg.gmm(p).labels.true.images)~=Np
            error(sprintf('Incompatible numbers of label images (pop-%d).',p));
        end
    end
end

cl  = cell(N,1);
dat = struct('dm',cl, 'Mat',cl, 'samp',[1 1 1], 'onam','', 'odir','',...
             'q',cl, 'v',cl, 'delta', zeros(1,K), 'psi',cl, 'model',cl, 'lab',cl, 'E',cl,'nvox',cl);
n   = 0;

% Process categorical data
if numel(cfg.cat)>=1
    Np = numel(cfg.cat{1});
    Nc = numel(cfg.cat);
    cl = cell(Nc,1);
    for np = 1:Np
        n = n + 1;
        for c=1:Nc
            cl{c} = cfg.cat{c}{np};
        end
        f         = nifti(char(cl));
        dm        = size(f(1).dat,[1 2 3]);
        dat(n).dm = dm;
        for c=2:Nc
            dmc = size(f(c).dat,[1 2 3]);
            if ~all(dmc==dm)
                error(sprintf('Incompatible image dimensions for tissue classes of subject %d (%dx%dx%d ~= %dx%dx%d)', np, dmc, dm));
            end
        end

        dat(n).Mat = f(1).mat;
        dat(n).q   = zeros(size(B,3),1);
        for c=2:Nc
            if sum((f(c).mat(:)-f(1).mat(:)).^2) > 1e-9 
                warning('Incompatible s-form matrices for subject %d (tissue classes)', np);
            end
        end

        [~,nam,~]   = fileparts(cl{1});
        dat(n).onam = sprintf('%d_%.5d_%s_%s', 0, np, nam, cfg.onam);
        dat(n).odir = sett.odir;
        dat(n).v    = fullfile(dat(n).odir,['v_' dat(n).onam '.nii']);
        dat(n).psi  = fullfile(dat(n).odir,['y_' dat(n).onam '.nii']);
        if isfinite(sett.del_settings)
            dat(n).delta = zeros(1,K);
        else
            dat(n).delta = [];
        end


        Kn = 0;
        for c=1:Nc
            Kn = Kn + size(f(1).dat,4);
        end
        if K>=0
            if Kn~=K
                error(sprintf('Incompatible numbers of categories for subject %d (%d~=%d).',np,Kn,K));
            end
        else
            K = Kn;
        end
        dat(n).model = struct('cat',struct('f',f, 'K',K));
        dat(n).lab   = [];
    end
end

% This is (for now) a hidden 'option' that enables to model
% multiple populations--with different number of observed
% channels--with the same GMM. The idea is that there are C
% channels in total, where each population (p) has C_p observed
% channels (it is required that at least one population has all C
% channels observed). Then, for populations that does not have all
% channels observed, their unobserved channels are given as images
% with all zeros/NaNs.
if isfield(cfg,'one_gmm_prior')
    one_gmm_prior = cfg.one_gmm_prior;
else
    one_gmm_prior = false;
end
num_priors = numel(cfg.gmm);
if one_gmm_prior
    num_priors = 1;
end

% Process scans (for gmm)
sett.gmm  = struct('pr',cell(num_priors,1),'hyperpriors',true, ....
                   'mg_ix', [], 'C',0, 'tol_gmm',[],'nit_gmm_miss',[],'nit_gmm',[],'nit_appear',[]);
for p=1:numel(cfg.gmm)

    ix_gmm = p;
    if one_gmm_prior
        ix_gmm = 1;
    end

    sett.gmm(ix_gmm).tol_gmm      = cfg.gmm(p).tol_gmm;
    sett.gmm(ix_gmm).nit_gmm_miss = cfg.gmm(p).nit_gmm_miss;
    sett.gmm(ix_gmm).nit_gmm      = cfg.gmm(p).nit_gmm;
    sett.gmm(ix_gmm).nit_appear   = cfg.gmm(p).nit_appear;

    % Multiple Gaussians per template class
    if isfield(cfg.gmm(p),'mg_ix') && ~isempty(cfg.gmm(p).mg_ix)
        sett.gmm(ix_gmm).mg_ix = cfg.gmm(p).mg_ix;
    else
        sett.gmm(ix_gmm).mg_ix = 1:(sett.K+1);
    end

    Nc = numel(cfg.gmm(p).chan);
    inu_reg  = zeros(Nc,1)+NaN;
    for c=1:Nc
        inu_reg(c)  = cfg.gmm(p).chan(c).inu.inu_reg;
    end
    sett.gmm(ix_gmm).inu_reg = inu_reg;

    cl = cell(Nc,1);
    C  = -1;
    if numel(cfg.gmm(p).chan)>=1
        Np = numel(cfg.gmm(p).chan(1).images);
        for np=1:Np
            n = n + 1;

            for c=1:Nc
                cl{c} = cfg.gmm(p).chan(c).images{np};
            end
            f          = nifti(char(cl));
            dm         = size(f(1).dat,[1 2 3]);
            dat(n).dm  = dm;
            dat(n).Mat = f(1).mat;
            dat(n).q   = zeros(size(B,3),1);

            for c=2:Nc
                dmc = size(f(c).dat,[1 2 3]);
                if ~all(dmc==dm)
                    error(sprintf('Incompatible image dimensions for subject %d in population %d (%dx%dx%d ~= %dx%dx%d)', np, p, dmc, dm));
                end
            end

            for c=2:Nc
                if sum((f(c).mat(:)-f(1).mat(:)).^2) > 1e-9
                    warning('Incompatible s-form matrices for subject %d in population %d', np, p);
                end
            end

            [~,nam,~]   = fileparts(cl{1});
            dat(n).onam = sprintf('%d_%.5d_%s_%s', p, np, nam, cfg.onam);
            dat(n).odir = sett.odir;
            dat(n).v    = fullfile(dat(n).odir,['v_' dat(n).onam '.nii']);
            dat(n).psi  = fullfile(dat(n).odir,['y_' dat(n).onam '.nii']);

            if isfinite(sett.del_settings)
                dat(n).delta = zeros(1,K);
            else
                dat(n).delta = [];
            end
            cf = zeros(Nc,1);
            for c=1:Nc
                cf(c) = size(f(1).dat,4);
            end
            Cn       = sum(cf);
            inu_co   = zeros(Cn,1);
            modality = zeros(Cn,1);
            ind      = 0;
            for c=1:Nc
                ind           = max(ind) + (1:cf(c));
                inu_co(ind)   = cfg.gmm(p).chan(c).inu.inu_co;
                modality(ind) = cfg.gmm(p).chan(c).modality;
            end

            % Set up INU fields
            T      = cell(1,Cn);
            vx     = sqrt(sum(dat(n).Mat(1:3,1:3).^2));
            for c=1:Cn
                if isnan(inu_co(c))
                    T{c} = zeros(0,0,0);
                elseif isinf(inu_co(c))
                    T{c} = zeros(1,1,1);
                else
                    d3   = ceil(2*vx.*dat(n).dm/inu_co(c));
                    T{c} = zeros(d3);
                end
            end

            if C>=0
                if Cn~=C
                    error(sprintf('Incompatible numbers of channels (pop=%d, n=%d).',p, np));
                end
            else
                C = Cn;
            end

            % Deal with labels
            if isfield(cfg.gmm(p).labels,'true')
               %lab = struct('f', nifti(cfg.gmm(p).labels.true.images{np}),...
               %             'cm_map', {cfg.gmm(p).labels.true.cm_map},...
               %             'w',       cfg.gmm(p).labels.true.w);
                lab = struct('f', nifti(cfg.gmm(p).labels.true.images{np}));
                dmc = size(lab.f(1).dat,[1 2 3]);
                if ~all(dmc==dm)
                    error(sprintf('Incompatible image dimensions for images of subject %d in population %d (%dx%dx%d ~= %dx%dx%d)', np, p, dmc, dm));
                end
                if sum((lab.f.mat(:)-f(1).mat(:)).^2) > 1e-9
                    warning('Incompatible s-form matrices for subject %d in population %d', np, p);
                end
               %if max(cellfun(@max,lab.cm_map)) > K+1 || min(cellfun(@min,lab.cm_map)) < 1
               %    error(sprintf('Poorly specified label mapping for population %d', p));
               %end
                Alpha = ones(256,K+1)*eps;
            else
                lab   = [];
                Alpha = [];
            end
            dat(n).lab   = lab;

            lb           = struct('sum', NaN, 'X', [], 'XB', [], 'Z', [], 'P', [], 'mu', [], 'A', [],'Alpha',[]);
            gmm          = struct('f',f, 'pop', ix_gmm, 'samp',[1 1 1],...
                                  'modality', modality, 'T',{T}, 'lb', lb,...
                                  'm',rand(Cn,K+1),'b',zeros(1,K+1)+1e-6,...
                                  'W',repmat(eye(Cn,Cn),[1 1 K+1]),'nu',zeros(1,K+1)+1e-6, 'gam',[],'Alpha',Alpha);
            dat(n).model = struct('gmm',gmm);
        end
    end

    % Load information from intensity priors file.
    % Note that such files would need to be hand-crafted.
    sett.gmm(ix_gmm).C     = C;
    sett.gmm(ix_gmm).pr    = {};
    sett.gmm(ix_gmm).hyperpriors = cfg.gmm(p).pr.hyperpriors;
    if ~isempty(cfg.gmm(p).pr.file) && ~isempty(cfg.gmm(p).pr.file{1})
        pr = load(cfg.gmm(p).pr.file{1});
        if isfield(pr,'mg_ix')
            if max(pr.mg_ix) ~= K+1
                error(sprintf('Incompatible K dimensions for intensity priors ("%s").',cfg.gmm(p).pr.file{1}));
            end
            sett.gmm(ix_gmm).mg_ix = pr.mg_ix;
        end
        if isfield(pr,'pr')
            if size(pr.pr{1},1) ~= C
                error(sprintf('Incompatible C dimensions for intensity priors ("%s").',cfg.gmm(p).pr.file{1}));
            end
            if size(pr.pr{1},2) ~= numel(sett.gmm(ix_gmm).mg_ix)
                error(sprintf('Incompatible total K dimensions for intensity priors ("%s").',cfg.gmm(p).pr.file{1}));
            end
            sett.gmm(ix_gmm).pr = pr.pr;
        end
    end
end

if isfield(sett.mu,'create')
    [sett.mu.Mmu,sett.mu.d] = specify_mean(dat,sett.mu.create.vx*[1 1 1]);
end
%==========================================================================

%==========================================================================
function [sett,dat] = random_init(sett,dat)
if isempty(sett.gmm), return; end

code = zeros(numel(dat),1);
for n=1:numel(dat)
    if isfield(dat(n).model,'gmm')
        code(n) = dat(n).model.gmm.pop;
    end
end

for p=1:numel(sett.gmm) % Loop over populations
    index = find(code==p); % Subjects in population p
    N     = numel(index);  % Number of subjects in population p
    C     = sett.gmm(p).C; % Number of channels for population p

    % Compute mean and variances for each subject in population p
    % and initialise INU parameters so subjects have the same mean
    mu_all = zeros(C,N); % Means
    vr_all = zeros(C,N); % Diagonal of covariance
    for n=1:N % Loop over subjects
        n1  = index(n);                   % Index of this subject
        gmm = dat(n1).model.gmm;          % GMM data for this subject
        dm  = dat(n1).dm;                 % Image dimensions
        f   = spm_mb_io('get_image',gmm); % Image data
        f   = reshape(f,prod(dm),C);      % Vectorise
        T   = gmm.T;                      % INU parameters
        mu  = zeros(C,1);                 % Mean
        vr  = zeros(C,1);                 % Diagonal of covariance
        for c=1:C                         % Loop over channels
            fc    = f(:,c);                        % Image for this channel
            fc    = fc(isfinite(fc));              % Ignore non-finite values
            if isempty(fc)
                T{c} = []; % No observations in channel => do not model bias field
            end
            mn    = min(fc);                      % Minimum needed for e.g. CT
            mu(c) = sum(fc)/size(f,1);            % Mean (assuming missing values are zero)
            fc    = fc(fc>((mu(c)-mn)/8+mn));     % Voxels above some threshold (c.f. spm_global.m)
            mu(c) = mean(fc);                     % Mean of voxels above the threshold
            vr(c) = var(fc);                      % Variance of voxels above the threshold
            m     = gmm.modality(c);              % Image modality
            if ~isempty(T{c}) && m~=2             % Should INU or global scaling be done?
                s           = 1000;               % Scale means to this value
                dc          = log(s)-log(mu(c));  % Log of scalefactor
                bbb         = spm_dctmtx(dm(1),1,1)*spm_dctmtx(dm(2),1,1)*spm_dctmtx(dm(3),1,1);
                T{c}(1,1,1) = dc/bbb;             % Adjust log-scalefactor to account for 3D DCT
                vr(c) = vr(c).*(s./mu(c)).^2;     % Adjust variance for rescaling
                mu(c) = s;                        % Set mean to s
            end
            vr_all(:,n) = vr; % Collect variances
            mu_all(:,n) = mu; % Collect means
        end
        dat(n1).model.gmm.T = T; % Assign INU parameters with new DC component
    end

    % Fill in entirely missing values with means
    for c=1:C
        msk  = isfinite(vr_all(c,:));
        vr_all(c,~msk) = sum(vr_all(c,msk),2)/sum(msk,2);
        msk  = isfinite(mu_all(c,:));
        mu_all(c,~msk) = sum(mu_all(c,msk),2)/sum(msk,2);
    end

    % Fix for when just a single subject, with one (or more) fully unobserved 
    % modalities
    mu_all(~isfinite(mu_all)) = 1000;
    vr_all(~isfinite(vr_all)) = mean(vr_all(isfinite(vr_all)));

    mg_ix  = sett.gmm(p).mg_ix; 
    K1     = numel(mg_ix); % Total number of Gaussians (some tissues may have more than one)
    gam    = ones(1,K1);
    Krange = unique(mg_ix);
    for k=Krange(:)'
        ind      = Krange==k;
        gam(ind) = 1/sum(ind);
    end

    if isempty(sett.gmm(p).pr)

        % If no priors specified, then generate some that are reasonably uninformative
        sett.gmm(p).pr    = cell(1,4);

        % Uninformative prior of mean - set to mean of image
        mu                = double(mean(mu_all,2));
        sett.gmm(p).pr{1} = repmat(mu,1,K1);
        sett.gmm(p).pr{2} = ones(1,K1)*1e-4;

        % Uninformative prior for variance - based on variance of image
        vr                = double(mean(vr_all,2));
        nu0               = C-1+1e-4;                % Minimally informative
        scale             = max(K1-1,1).^(2/C);      % Crude heuristic
        W0                = diag(1./vr)*(scale/nu0);
        sett.gmm(p).pr{3} = repmat(W0,[1 1 K1]);
        sett.gmm(p).pr{4} = ones(1,K1)*nu0;

        % Uninformative priors for gam - unused
        sett.gmm(p).pr{5} = ones(1,K1)*eps;

        % Uninformative priors for Alpha - used for categorical mixture
        sett.gmm(p).pr{6} = ones(256,K1)*eps; 

        % Random mean intensities, roughly sorted. Used to break symmetry.
        rng('default'); rng(1); % Want some reproducibility
        mu                = bsxfun(@plus,0.01*diag(sqrt(vr)*(1-1/scale))*randn(C,K1), mu);
        d                 = sum(diag(sqrt(vr*(1-1/K1)))\mu,1);           % Heuristic measure of how positive
        [~,o]             = sort(-d); % Order the means, most positive first
        mu                = mu(:,o);

        % Assign the same GMM starting estimates for all subjects
        for n=1:N
            n1                    = index(n);
            dat(n1).model.gmm.m   = mu;                % Random means (break symmetry)
           %dat(n1).model.gmm.m   = bsxfun(@plus,0.01*diag(sqrt(vr)*(1-1/scale))*randn(C,K1), mu);
            dat(n1).model.gmm.b   = sett.gmm(p).pr{2};
            dat(n1).model.gmm.W   = sett.gmm(p).pr{3};
            dat(n1).model.gmm.nu  = sett.gmm(p).pr{4};
            dat(n1).model.gmm.gam = gam;
        end
    else
        for n=1:N
            n1   = index(n);

%           % Initial distribution for mean
%           m    = sett.gmm(p).pr{1};      % Use prior mean
%           b    = ones(1,K1)*1e-3;        % Uninformative
%
%           % Initial distribution for precision
%           nu0  = size(m,1)-1+1e-3;
%           nu   = ones(1,K1)*nu0;
%
%           vr   = double(mean(vr_all,2));
%           scal = max(K1-1,1).^(2/C);     % Crude heuristic
%           W    = diag(1./vr)*(scal/nu0); % Low precision
%           W    = repmat(W,[1 1 K1]);

            % directly using the learned priors seems to work better
            % (for CT at least)
            m = sett.gmm(p).pr{1};
            b = sett.gmm(p).pr{2};
            W = sett.gmm(p).pr{3};
            nu = sett.gmm(p).pr{4};

            dat(n1).model.gmm.m   = m;
            dat(n1).model.gmm.b   = b;
            dat(n1).model.gmm.W   = W;
            dat(n1).model.gmm.nu  = nu;
            dat(n1).model.gmm.gam = gam;
            if isa(dat(n1).lab,'struct')
                dat(n1).model.gmm.Alpha = ones(256,K1)*eps;
            else
                dat(n1).model.gmm.Alpha = [];
            end
        end
    end
end
%==========================================================================

%==========================================================================
function [Mat,dm] = specify_mean(dat,vx)
% Specify the voxel-to-world matrix (Mat) and dimensions (dm) of the template
N  = numel(dat);     % Number of subjects
if N==0
    Mat = eye(4);
    dm  = [0 0 0];
    return
end
mn =  Inf*ones(3,N); % Minimum coordinates
mx = -Inf*ones(3,N); % Maximum coordinates
for n=1:N
    dm      = [dat(n).dm 1 1];
    corners = [
        1 dm(1)    1  dm(1)   1  dm(1)    1  dm(1)
        1    1  dm(2) dm(2)   1     1  dm(2) dm(2)
        1    1     1     1 dm(3) dm(3) dm(3) dm(3)
        1    1     1     1    1     1     1     1]; % Voxel indices at corners
    M       = dat(n).Mat;
    pos     = M(1:3,:)*corners; % mm coordinates of corners
    mx(:,n) = max(pos,[],2);    % Maximum mm coordinates for subject n
    mn(:,n) = min(pos,[],2);    % Minimum mm coordinates for subject n
end

% Final bounding box should enclose 95% of the scans
ind = min(max(round(N*0.95),1),N); % 95th percentile
mx  = sort(mx,2,'ascend');
mn  = sort(mn,2,'descend');
mx  = mx(:,ind);
mn  = mn(:,ind);
bb  = [mn'; mx'];

vx  = [-1 1 1] .* abs(vx);
%mn  = vx .* min(bb ./ vx);        % "first" voxel's mm coordinates
mn = bsxfun(@times,vx,min(bsxfun(@rdivide,bb,vx)));
%mx  = vx .* round(max(bb ./ vx)); % "last voxel's mm coords
mx = bsxfun(@times,vx,max(bsxfun(@rdivide,bb,vx)));
Mat = spm_matrix([mn 0 0 0 vx]) * spm_matrix([-1 -1 -1]);
dm  = Mat \ [mx 1]';
dm  = round(dm(1:3)');
%==========================================================================
