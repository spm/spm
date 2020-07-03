function varargout = spm_mb_appearance(action,varargin) % Appearance model
%
% FORMAT chan       = spm_mb_appearance('inu_basis',T,df,Mat,reg,samp)
% FORMAT [inu,ll]   = spm_mb_appearance('inu_field',T,chan,d,varargin)
% FORMAT labels     = spm_mb_appearance('get_labels',dat,sett,do_samp)
% FORMAT z          = spm_mb_appearance('responsibility',m,b,V,n,f,mu,msk_chn)
% FORMAT [z,dat]    = spm_mb_appearance('update',dat,mu0,sett)
% FORMAT dat        = spm_mb_appearance('update_prior',dat,sett)
%__________________________________________________________________________
% Copyright (C) 2019-2020 Wellcome Centre for Human Neuroimaging

% $Id: spm_mb_appearance.m 7885 2020-07-03 14:10:31Z mikael $


switch action
    case 'inu_basis'
        [varargout{1:nargout}] = inu_basis(varargin{:});
    case 'inu_field'
        [varargout{1:nargout}] = inu_field(varargin{:});
    case 'inu_recon'
        [varargout{1:nargout}] = inu_recon(varargin{:});
    case 'get_labels'
        [varargout{1:nargout}] = get_labels(varargin{:});
    case 'responsibility'
        [varargout{1:nargout}] = responsibility(varargin{:});
    case 'update'
        [varargout{1:nargout}] = update(varargin{:});
    case 'update_prior'
        [varargout{1:nargout}] = update_prior(varargin{:});
    case 'update_shared_prior'
        [varargout{1:nargout}] = update_shared_prior(varargin{:});
    otherwise
        error('Unknown function %s.', action)
end
%==========================================================================

%==========================================================================
function [inu,ll] = inu_field(T,chan)
d  = [size(chan(1).B1,1) size(chan(1).B2,1) size(chan(1).B3,1)];
nz = d(3);
C  = numel(T);
cr = 1:C;

% Compute full bias field (for all channels)
inu = zeros([d C],'single');
ll  = zeros(1,C);

for c=cr
    t = double(T{c});
    ll(c) = -0.5*t(:)'*chan(c).ICO*t(:);
    for z=1:nz
        inu_c        = inu_transform(chan(c).B1,chan(c).B2,chan(c).B3(z,:),T{c});
        inu(:,:,z,c) = single(exp(inu_c));
    end
end
%==========================================================================

%==========================================================================
function t = inu_transform(B1,B2,B3,T)
% Create an image-space log bias field from its basis function encoding.
%
% FORMAT t = inu_transform(B1,B2,B3,T)
% B1 - x-dim DCT basis [nx kx]
% B2 - y-dim DCT basis [ny ky]
% B3 - z-dim DCT basis [nz kz]
% T  - DCT encoding of the log bias field [kx ky kz]
% t  - Reconstructed log bias field [nx ny nz]
if ~isempty(T)
    d2 = [size(T) 1];
    t1 = reshape(reshape(T, d2(1)*d2(2),d2(3))*B3', d2(1), d2(2));
    t  = B1*t1*B2';
else
    t  = zeros(size(B1,1),size(B2,1));
end
%==========================================================================

%==========================================================================
function chan = inu_basis(T,df,Mat,reg,samp)
if nargin<5, samp = 0; end
if nargin<4, reg  = ones(1,numel(T)); end
cl   = cell(1, numel(T));
chan = struct('ICO', cl, 'B1',cl, 'B2',cl, 'B3',cl);
ind  = sample_ind(df,samp);
vx   = sqrt(sum(Mat(1:3,1:3).^2,1));
for c=1:numel(T)
    d3 = [size(T{c}) 1];
    d3 = d3(1:3);

    % GAUSSIAN REGULARISATION for bias correction
    chan(c).ICO = reg(c)*spm_bias_lib('regulariser','bending',df,d3,vx);

    % Basis functions for bias correction
    chan(c).B1 = spm_dctmtx(df(1),d3(1),ind{1});
    chan(c).B2 = spm_dctmtx(df(2),d3(2),ind{2});
    chan(c).B3 = spm_dctmtx(df(3),d3(3),ind{3});
end
%==========================================================================

%==========================================================================
function [tr,D] = inu_reg_suffstat(T,Sig,df,Mat)
vx   = sqrt(sum(Mat(1:3,1:3).^2,1));
C    = numel(T);
tr   = zeros(1,C);
D    = zeros(1,C);
for c=1:numel(T)
    d3   = size(T{c});
    D(c) = prod(d3); % Should really subtract 1 because of the DC.
    if prod(d3)>0
        ICO = spm_bias_lib('regulariser','bending',df,d3,vx);
        t   = T{c}(:);
        S   = Sig{c};
        if size(S,2)==1
            % Diagonal approximation
            tr(c) = sum(S.*t.^2) + sum(S.*diag(IC0));
        else
            tr(c) = t'*ICO*t + sum(sum(S.*ICO));
        end
    end
end
%==========================================================================

%==========================================================================
function labels = get_labels(dat,sett,samp)
if nargin < 3, samp = [1 1 1]; end
lab = dat.model.gmm.lab;
if isempty(lab), labels = 0; return; end

cm_map     = lab.cm_map; % cell array that defines the confusion matrix (cm)

% Load labels
labels = spm_mb_io('get_data',lab.f);
if any(samp~=1)
    labels = subsample(labels,samp);
end

% Use labels2use to keep only labels of interest
labels = round(labels(:));
labels(labels<1 | labels>numel(cm_map)) = numel(cm_map);

% Get confusion matrix that maps from label value to probability value
cm = get_label_conf_matrix(cm_map,lab.w,sett.K);

% Build NxK1 label image using confusion matrix
labels = cm(labels,:);
%==========================================================================

%==========================================================================
function cm = get_label_conf_matrix(cm_map,w,K)
% FORMAT CM = get_label_cm(cm_map,opt)
% cm_map - Defines the confusion matrix
% sett   - Options structure
% cm     - confusion matrix
%
% Build Rater confusion matrix for one subject.
% This matrix maps template classes to manually segmented classes.
% Manual labels often do not follow the same convention as the Template,
% and not all regions may be labelled. Therefore, a manual label may
% correspond to several Template classes and, conversely, one Template
% class may correspond to several manual labels.
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Here, we assume that all subjects from the same population (e.g.,
% a publicily available dataset) have the same labelling protocole and
% confusion matrix.
% We allow the rater's sensitivity to change every few acquistion. We would
% typically start with a high sensitivity, to weight the labels strongly,
% and then decrease this value to allow the model to correct the rater's
% mistakes (especially near boundaries).

% Parse function settings
w  = min(max(w,1e-7),1-1e-7);
K1 = K + 1;
L  = numel(cm_map); % Number of labels
cm = zeros([L K1],'single'); % Allocate confusion matrix
for l=1:L % Loop over labels
    ix            = false(1,K1);
    ix(cm_map{l}) = true;
    cm(l, ix)     = log(w/sum( ix));
    cm(l,~ix)     = log((1-w)/sum(~ix));
end
%==========================================================================

%==========================================================================
function [z,lb] = responsibility(m,b,V,n,mf,vf,mu,msk_chn)
% Compute responsibilities.
%
% FORMAT z = responsibility(m,b,V,n,f,vf,mu,L,code)
% m       - GMM Means
% b       - GMM Mean d.f.
% V       - GMM Scale matrices
% n       - GMM Scale d.f.
% mf      - Expectation of INU-corrected observed image in matrix form [nbvox nbchannel]
% vf      - Variance of INU-corrected observed image in matrix form [nbvox nbchannel]
% mu      - Deformed and exponentiated template
% msk_chn - Mask of observed channels per code
% z       - Image of responsibilities [nbvox K]

const  = spm_gmm_lib('Normalisation', {m,b}, {V,n}, msk_chn);
if ~isempty(vf)
    z      = spm_gmm_lib('Marginal', mf, {m,V,n}, const, msk_chn, vf);
else
    z      = spm_gmm_lib('Marginal', mf, {m,V,n}, const, msk_chn);
end
[z,lb] = spm_gmm_lib('Responsibility', z, mu);
%==========================================================================

%==========================================================================
function [z,lb] = responsibility_t(m,b,V,n,mf,vf,mu,msk_chn)
% Compute responsibilities.
%
% FORMAT z = responsibility(m,b,V,n,f,vf,mu,L,code)
% m       - GMM Means
% b       - GMM Mean d.f.
% V       - GMM Scale matrices
% n       - GMM Scale d.f.
% mf      - Expectation of INU-corrected observed image in matrix form [nbvox nbchannel]
% vf      - Variance of INU-corrected observed image in matrix form [nbvox nbchannel]
% mu      - Deformed and exponentiated template
% msk_chn - Mask of observed channels per code
% z       - Image of responsibilities [nbvox K]

if ~isempty(vf)
    z      = spm_gmm_lib('Marginal_t', mf, {m,b,V,n}, msk_chn,vf);
else
    z      = spm_gmm_lib('Marginal_t', mf, {m,b,V,n}, msk_chn);
end
[z,lb] = spm_gmm_lib('Responsibility', z, mu);
%==========================================================================

%==========================================================================
function [dat,Z] = update(dat,mu0,sett)
% Update appearance model for a single subject (GMM & bias field)
%
% FORMAT [dat,Z] = update(dat,mu0,sett)
% dat - Structure holding data for a single subject
% mu0 - Log template
% sett - Structure of settings

% Parse function settings
gmm          = dat.model.gmm;
mg_ix        = sett.gmm(gmm.pop).mg_ix;
nit_gmm      = sett.gmm(gmm.pop).nit_gmm;
nit_gmm_miss = sett.gmm(gmm.pop).nit_gmm_miss;
nit_appear   = sett.gmm(gmm.pop).nit_appear;
tol_gmm      = sett.gmm(gmm.pop).tol_gmm;
C            = sett.gmm(gmm.pop).C;
inu_reg      = sett.gmm(gmm.pop).inu_reg;
Kmg          = numel(mg_ix);
Mat          = dat.Mat;
df           = dat.dm;
ds           = [size(mu0) 1 1];
ds           = ds(1:3);
K1           = sett.K + 1;

% Get image data
f0    = spm_mb_io('get_image',gmm);
samp  = dat.samp;
samp2 = dat.samp2;
f0    = subsample(f0,dat.samp);
% Labels are of size [df0 K+1], so cannot be added to template at this stage
lab0  = get_labels(dat,sett,samp);
if lab0 ~= 0
    lab0  = vec2vol(lab0,[size(f0,1) size(f0,2) size(f0,3) K1]);
end

% Intensity priors
pr   = sett.gmm(gmm.pop).pr;
gmm  = fix_scaling(gmm,pr,df);

% GMM posterior
m    = gmm.m;
b    = gmm.b;
V    = gmm.V;
n    = gmm.n;
mg_w = gmm.mg_w;

% % Broaden the variance to make less informative.
% % Might sometimes help escape local optima. Broadened
% % more when sample density is lowest. Needs further
% % testing.
% scal = prod(samp.*samp2).^(-1/4);
% b    = b*scal;
% n    = (n-C)*scal+C;

% If template is missing, set corresponding voxels of
% f to missing too.
% Note that if the template doesn't cover the entire field
% of view, this can cause the log-likelihood to behave oddly.
msk_vx = any(isnan(mu0),4);
for c=1:size(f0,4)
    fc = f0(:,:,:,c);
    fc(msk_vx) = NaN;
    f0(:,:,:,c) = fc;
end
clear msk_vx

% Compute Gaussian parameters on a subset of the voxels.
% Note strange behaviour for interleaved images where there
% are systematic differences between odd and even slices.
% If the subset contains only odd slices, then the overall
% log-likelihood can increase.
[f,d] = subsample( f0,samp2);
mu    = subsample(mu0,samp2);
lab   = vol2vec(subsample(lab0,samp2));

% Template
mu    = vol2vec(spm_mb_shape('template_k1',mu)); % Make K + 1 template
mu    = bsxfun(@plus,mu,lab); % Add labels to template
clear lab
mu    = mu(:,mg_ix); % Expand, if using multiple Gaussians per tissue

% Bias field related
T     = gmm.T;
if isfield(gmm,'Sig')
    Sig    = gmm.Sig;
else
    Sig    = cell(size(T));
end

do_inu = ~cellfun(@isempty,T);
if any(do_inu)
    chan          = inu_basis(T,df,Mat,inu_reg,samp.*samp2);
    [llinu,mf,vf] = inu_recon(f,chan,T,Sig);
    lxb           = sum(llinu(:),'double');
else
    mf            = f;
    vf            = zeros(size(f),'single');
    lxb           = 0;
end

% Format for spm_gmm
[mf,code_image,msk_chn] = spm_gmm_lib('obs2cell', vol2vec(mf));
mu                      = spm_gmm_lib('obs2cell', mu, code_image, false);
vf                      = spm_gmm_lib('obs2cell', vol2vec(vf), code_image, true);
code_list               = unique(code_image);
code_list               = code_list(code_list ~= 0);
nvox                    = sum(code_image(:)>0);
lbs = -Inf;

for it_appear=1:nit_appear

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update GMM and get responsibilities (Z)
    %------------------------------------------------------------
    lbso = lbs;
    [Z,mog,~,lb,mg_w] = spm_gmm_lib('loop',mf,...
                                    {{m,b},{V,n}},{'LogProp', mu}, ...
                                    'ObsUncertainty', vf, ...
                                    'GaussPrior',     pr, ...
                                    'Missing',        msk_chn, ...
                                    'IterMax',        nit_gmm, ...
                                    'Tolerance',      tol_gmm*nvox, ...
                                    'SubIterMax',     nit_gmm_miss, ...
                                    'SubTolerance',   tol_gmm*nvox*0.1, ...
                                    'Verbose',        [0 0], ...
                                    'MultGaussPi',    {mg_ix,mg_w});
    m      = mog.mu;
    b      = mog.b;
    V      = mog.V;
    n      = mog.n;
    lx     = lb.X;
    lbs    = sum(lb.mu,'double')+sum(lb.A,'double')+lx+lxb;
    %fprintf('%g ', lbs/nvox);
    if (it_appear==nit_appear) || (lbs-lbso < tol_gmm*nvox) 
        % Finished
        break
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update bias field parameters
    % This computes the derivatives of the negative logarithm of the
    % joint probability distribution.
    % The uncertainty with which the INU is estimated could be
    % accounted for by extending the variational Bayes framework
    % to account for the estimation of the INU. This could use:
    %   E[f \exp(b)] = f \exp(\mu_b+\sigma^2_b/2)
    % Var[f \exp(b)] = f^2 \exp(2 \mu_b + \sigma^2_b) (\exp(\sigma^2_b)-1)
    %------------------------------------------------------------
    if any(do_inu)

        % Update bias field parameters for each channel separately
        for c=1:C % Loop over channels
            if isempty(T{c}), continue; end

            % Compute mode, rather than expectations for the parameters
            % to update.
            [~,mf_c,vf_c] = inu_recon(f(:,:,:,c),chan(c),T(c));
            mf = insert2cell(mf,code_image,vol2vec(mf_c),c);
            vf = insert2cell(vf,code_image,vol2vec(vf_c),c);

            % Compute gradient and Hessian (in image space)
            gr_im = zeros(d(1:3),'single');
            H_im  = zeros(d(1:3),'single');
            for l=1:size(msk_chn,1) % loop over combinations of missing voxels

                % Get mask of missing modalities (with this particular code)
                ixo = msk_chn(l,:);         % Observed channels
                ixm = ~ixo;                 % Missing channels
                nm  = sum(ixm);             % Number of missing channels
                if ~ixo(c), continue; end

                % Convert channel indices to observed indices
                ixc = 1:C; % mapped_c
                ixc = ixc(ixo);
                ixc = find(ixc == c);

                go = 0; % Gradient accumulated accross clusters
                Ho = 0; % Hessian accumulated accross clusters
                for k=1:Kmg
                    % Compute expected precision (see GMM + missing data)
                    Voo = V(ixo,ixo,k);
                    Vom = V(ixo,ixm,k);
                    Vmm = V(ixm,ixm,k);
                    Vmo = V(ixm,ixo,k);
                    Ao  = Voo - Vom*(Vmm\Vmo);
                    Ao  = (n(k) - nm) * Ao;
                    mo  = m(ixo,k);

                    % Compute statistics
                    gk = bsxfun(@minus, mf{l}, mo.') * Ao(ixc,:).';
                    Hk = Ao(ixc,ixc);

                    gk = bsxfun(@times, gk, Z{l}(:,k));
                    Hk = bsxfun(@times, Hk, Z{l}(:,k));

                    % Accumulate across clusters
                    go = go + gk;
                    Ho = Ho + Hk;
                end

                % Multiply with bias corrected value (chain rule)
                omf = mf{l}(:,ixc);
                go  = go .*  omf;
                Ho  = Ho .* (omf.^2);
                clear omf

                % Add terms related to the normalisation (log(b))
                go = go - 1;
                Ho = Ho + 1; % Comes from subs(H,g,0)

                % Accumulate across missing codes
                ixvx        = (code_image == code_list(l));
                gr_im(ixvx) = gr_im(ixvx) + go;
                H_im(ixvx)  = H_im(ixvx)  + Ho;
                clear ixvx
            end
            clear Z

            % Actual Hessian is greater than the expected Hessian at some voxels,
            % so take the maximum of the expected and actual.
            H_im(gr_im<0) = H_im(gr_im<0) - gr_im(gr_im<0);

            % Compute gradient and Hessian
            d3 = numel(T{c}); % Number of DCT parameters
            H  = zeros(d3,d3);
            gr = zeros(d3,1);

            B1 = double(chan(c).B1);
            B2 = double(chan(c).B2);
            B3 = double(chan(c).B3);

            for z=1:d(3)
                b3 = B3(z,:)';
                gr = gr + kron(b3,spm_krutil(double(gr_im(:,:,z)),B1,B2,0));
                H  = H  + kron(b3*b3',spm_krutil(double(H_im(:,:,z)),B1,B2,1));
            end
            clear gr_im H_im b3

            % Gauss-Newton update of bias field parameters
            S      = inv(H + chan(c).ICO);
            Sig{c} = diag(S);
            T{c}   = T{c} - reshape(S*(gr + chan(c).ICO*T{c}(:)),size(T{c}));
            clear H gr

            % Compute new expectations (only for channel c)
            [llinu(:,c),mf_c,vf_c] = inu_recon(f(:,:,:,c),chan(c),T(c),Sig(c));
            mf     = insert2cell(mf,code_image,vol2vec(mf_c),c);
            vf     = insert2cell(vf,code_image,vol2vec(vf_c),c);
            lxb    = sum(llinu(:),'double');
            Z      = responsibility(m,b,V,n,mf,vf,reweight_mu(mu,log(mg_w)),msk_chn);
        end
    end
end
clear f mu

% Update dat
lbs      = lx+lxb+sum(lb.mu,'double')+sum(lb.A,'double');
gmm.T    = T;
gmm.Sig  = Sig;
gmm.m    = m;
gmm.b    = b;
gmm.V    = V;
gmm.n    = n;
gmm.lb   = lb;
gmm.mg_w = mg_w;

if nargout > 1

    % Compute full size resps
    if any(samp2~=1)
        msk = subsample_mask(df,samp2);

        % Compute full-sized responsibilities on original data
        if any(do_inu) % Bias correct
            chan          = inu_basis(T,df,Mat,inu_reg,samp);
            [llinu,mf,vf] = inu_recon(f0,chan,T,Sig);
            lxb = sum(llinu(:));
        else
            mf  = f0;
            lxb = 0;
        end

        % Template
        mu0 = spm_mb_shape('template_k1',mu0);
        mu0 = bsxfun(@plus,mu0,lab0); % Add labels to template
        mu0 = bsxfun(@plus, mu0(:,:,:,mg_ix), reshape(log(mg_w),[1 1 1 numel(mg_w)]));

       %[mf1,code_image,msk_chn] = spm_gmm_lib('obs2cell', mask(msk,mf));
       %mu1     = spm_gmm_lib('obs2cell', mask(msk,mu0), code_image, false);
       %vf1     = spm_gmm_lib('obs2cell', mask(msk, vf), code_image, true);
       %[Z1,lx] = responsibility(m,b,V,n,mf1,vf1,mu1,msk_chn);
        Z       = unmask1(msk,collapse_Z(spm_gmm_lib('cell2obs', Z, code_image, msk_chn),mg_ix));

        % Compute other responsibilities from a mixture of Student's t distributions.
        % See Eqns. 10.78-10.82 & B.68-B.72 in Bishop's PRML book.
        % In practice, it only improves probabilities by a tiny amount.
        msk      = ~msk;
        [mf1,code_image,msk_chn] = spm_gmm_lib('obs2cell', mask(msk,mf));
        mu1      = spm_gmm_lib('obs2cell', mask(msk,mu0), code_image, false);
        vf1      = spm_gmm_lib('obs2cell', mask(msk, vf), code_image, true);
        [Z2,lx2] = responsibility_t(m,b,V,n,mf1,vf1,mu1,msk_chn);
        Z        = unmask1(msk,collapse_Z(spm_gmm_lib('cell2obs', Z2, code_image, msk_chn),mg_ix),Z);
        nvox     = nvox + sum(code_image(:)>0);
        lbs      = lx+lx2+lxb+lb.mu+lb.A;
        clear mu0
    else
        Z = vec2vol(collapse_Z(spm_gmm_lib('cell2obs', Z, code_image, msk_chn),mg_ix),ds);
    end
end
dat.E(1)      = -lbs; %*scl_samp;
dat.nvox      = nvox;
dat.model.gmm = gmm;
%fprintf(' [%g] ', dat.E(1)/dat.nvox);
%==========================================================================

%==========================================================================
function Z = collapse_Z(Z,mg_ix)
% If using multiple Gaussians per tissue, collapse so that Z is of
% size K.
% Note: This assumes that mg_ix is sorted into increasing values.
K = max(mg_ix)-1;
if size(Z,2) > K
    for k=1:K
        Z(:,k) = sum(Z(:,mg_ix==k),2);
    end
end
Z = Z(:,1:K);
%==========================================================================

%==========================================================================
function X2d = vol2vec(X4d)
d   = [size(X4d) 1 1 1];
X2d = reshape(X4d,[prod(d(1:3)) d(4)]);
%==========================================================================

%==========================================================================
function X4d = vec2vol(X2d,dm)
dm = [dm(:)' 1 1];
if size(X2d,1)~=prod(dm(1:3))
    error('Incompatible dimensions.');
end
X4d = reshape(X2d,[dm(1:3) size(X2d,2)]);
%==========================================================================

%==========================================================================
function gmm = fix_scaling(gmm,pr,df)
% This function appears to make things worse for some reason. This will need
% more thought to understand the reason why.

% Adjust the DC part of INU to best match the GMM posteriors to the priors.
msk    = ~cellfun(@isempty,gmm.T(:));    % Only rescale when required
msk    = msk & any(diff(pr{1},[],2),2);  % Only rescale for informative priors
if ~any(msk), return; end                % Return if nothing to do

po     = {gmm.m, gmm.b, gmm.V, gmm.n};
r      = zeros(size(pr{1},1),1);
r(msk) = get_scaling(pr,po,msk);
s      = exp(r);

% Adjust distributions of mean and precision matrices
gmm.m  = diag(1./s)*gmm.m;
for k=1:size(gmm.V,3)
    gmm.V(:,:,k) = diag(s)*gmm.V(:,:,k)*diag(s);
end

% Make corresponding change to bias field (accounting for DCT scaling)
for m=find(msk(:)')
    gmm.T{m}(1,1,1) = gmm.T{m}(1,1,1) - r(m)*sqrt(prod(df));
end
%==========================================================================

%==========================================================================
function r = get_scaling(pr,po,msk)
% Determine log of rescaling factor that best matches the priors
% with the posteriors. This was based on Eq. 10.74 of Bishop's
% PRML book.
% Substitute m0 in 10.74 for diag(exp(r))*mu0, and V0 for
% diag(exp(r))\V0/diag(exp(r)). Then differentiate w.r.t. r
% to obtain gradients and Hessians for a Newton optimisation.

[mu0,b0,V0,nu0] = deal(pr{:});
[mu ,~ ,V ,nu ] = deal(po{:});
[D,K]           = size(mu0);
if nargin<3, msk = true(D,1); end

% Ad hoc fix to reduce the chance of one class dominating the scaling.
% Perhaps some form of hyper-priors might be a more elegant solution.
b0   = min(b0,exp(mean(log(b0)))*10); 

Alph = 0;
beta = 0;
gamm = 0;
for k=1:K
    Alph = Alph + nu(k)*(inv(V0(:,:,k)).*V(:,:,k)' + ...
                  b0(k)*diag(mu0(:,k))'*V(:,:,k)*diag(mu0(:,k)));
    beta = beta + b0(k)*nu(k)*diag(mu0(:,k))'*V(:,:,k)*mu(:,k);
    gamm = gamm + nu0(k);
end
Alph = Alph(msk,msk);
beta = beta(msk);
r    = zeros(sum(msk),1);
for it=1:100
    s  = exp(r);
    g  = s.*(Alph*s - beta -  gamm./(s   +eps));
    H  = (s*s').*(Alph + diag(gamm./(s.^2+eps)));
    H  = H + eye(size(H))*1e-6*max(diag(H));
    dr = H\g;
    r  = r - dr;
    if sqrt(mean(dr.^2)) < 1e-9, break; end
end
%==========================================================================

%==========================================================================
function sett = update_prior(dat, sett)
if isempty(sett.gmm), return; end

% Get population indices
code = zeros(numel(dat),1);
for n=1:numel(dat)
    if isfield(dat(n).model,'gmm')
        code(n) = dat(n).model.gmm.pop;
    end
end

for p=1:numel(sett.gmm) % Loop over populations
    if iscell(sett.gmm(p).hyperpriors)
        index = find(code==p);
        N     = numel(index);
        pr    = sett.gmm(p).pr;

        % Get all posteriors
        po      = cell(1,N);
        for n=1:N
            n1          = index(n);
            po{n}{1}{1} = dat(n1).model.gmm.m;
            po{n}{1}{2} = dat(n1).model.gmm.b;
            po{n}{2}{1} = dat(n1).model.gmm.V;
            po{n}{2}{2} = dat(n1).model.gmm.n;
        end

        % Update prior
        hp = sett.gmm(p).hyperpriors;
        sett.gmm(p).pr = spm_gmm_lib('updatehyperpars',po,pr,hp{:});

        %% Update INU regularisation. Disabled because it under-regularises
        %ss_inu0 = zeros(1,size(pr{1},1));
        %ss_inu1 = zeros(1,size(pr{1},1));
        %for n=1:N
        %    [tr,D] = inu_reg_suffstat(dat(n1).model.gmm.T,dat(n1).model.gmm.Sigma,dat(n1).df,dat(n1).Mat);
        %    ss_inu1 = ss_inu1 + tr;
        %    ss_inu0 = ss_inu0 + D;
        %end
        %sett.gmm(p).inu_reg = ss_inu1./ss_inu2;

    end
end
%==========================================================================

%==========================================================================
function sett = update_shared_prior(dat, sett)
if isempty(sett.gmm), return; end

% Get population indices
code = zeros(numel(dat),1);
for n=1:numel(dat)
    if isfield(dat(n).model,'gmm')
        code(n) = dat(n).model.gmm.pop;
    end
end

for p=1:numel(sett.gmm) % Loop over populations
    if iscell(sett.gmm(p).hyperpriors)
        index = find(code==p);
        N     = numel(index);
        pr    = sett.gmm(p).pr;
        K     = size(pr{1},2);
        pr{1} = mean(pr{1},2);
        pr{2} = mean(pr{2},2);
        pr{3} = mean(pr{3},3);
        pr{4} = mean(pr{4},2);

        % Get all posteriors
        po    = cell(1,N*K);
        for n=1:N
            n1  = index(n);
            gmm = dat(n1).model.gmm;
            for k=1:K
                nk     = (n-1)*K+k;
                po{nk} = {{gmm.m(:,k),gmm.b(:,k)},{diag(diag(gmm.V(:,:,k))),gmm.n(:,k)}};
            end
        end

        % Update prior
        hp    = sett.gmm(p).hyperpriors;
        pr    = spm_gmm_lib('updatehyperpars',po,pr,hp{:});
        pr{1} = repmat(pr{1},1,K);
        pr{2} = repmat(pr{2},1,K);
        pr{3} = repmat(pr{3},1,1,K);
        pr{4} = repmat(pr{4},1,K);
        sett.gmm(p).pr = pr;
    end
end
%==========================================================================

%==========================================================================
%
% Utility functions
%
%==========================================================================

%==========================================================================
function mu = reweight_mu(mu,logmg_w)
if sum(logmg_w) == 0, return; end
for i=1:numel(mu)
   %mu{i} = mu{i} + logmg_w;
    mu{i} = bsxfun(@plus, mu{i}, logmg_w);
end
%==========================================================================

%==========================================================================
function [of,d,w] = subsample(f,samp)
% Subsample a multichannel volume.
%
% FORMAT [of,d,scl_samp] = subsample(f,samp);
% f    - Original volume
% samp - Sampling distances in voxels
% of   - Resampled volume
% d    - Output dimensions
% w    - Proportion of sampled voxels

if all(samp==1)
    of = f;
    d  = [size(f,1) size(f,2) size(f,3)];
    w  = 1;
else
    % Input image properties
    df   = [size(f) 1];
    df   = df(1:3);

    ind  = sample_ind(df,samp);
    d    = cellfun(@length,ind);  % New dimensions
    of   = f(ind{:},:,:);

    % For weighting data parts of lowerbound with factor based on amount of
    % downsampling
    w    = prod(df)/prod(d);
end
%==========================================================================

%==========================================================================
function msk = subsample_mask(d,samp)
% Generate a binary mask of which voxels were sampled
msk = false(d);
ind = sample_ind(d,samp);
msk(ind{:}) = true;
%==========================================================================

%==========================================================================
function varargout = mask(msk,varargin)
nt  = sum(msk(:));
dm  = size(msk);
if dm(2)==1, dm = dm(1); end % MATLAB fudge because ndim(msk)>=2
[ind{1:numel(dm)}] = deal(':');
varargout = cell(1,nargin-1);
for n=1:nargin-1
    vol_in  = varargin{n};
    d       = [size(vol_in) 1];
    dt      = d((numel(dm)+1):end);
    vec_out = zeros([nt dt],'like',vol_in);
    for n1=1:prod(dt)
        tmp = vol_in(ind{:},n1);
        vec_out(:,n1) = tmp(msk(:));
    end
    varargout{n} = vec_out;
end
%==========================================================================

%==========================================================================
function vol_out = unmask1(msk,vec_in,vol_out)
dm  = size(msk);
if dm(2)==1, dm = dm(1); end % MATLAB fudge because ndim(msk)>=2
[ind{1:numel(dm)}] = deal(':');
d       = size(vec_in);
dt      = d(2:end);
if nargin<3
    vol_out = zeros([dm dt],'like',vec_in);
end
for n1=1:prod(dt)
    tmp         = zeros([dm 1],'like',vec_in);
    tmp(msk(:)) = vec_in(:,n1);
    vol_out(ind{:},n1) = tmp;
end
%==========================================================================

%==========================================================================
function ind = sample_ind(df,samp)
sk   = max([1 1 1],samp);
ind  = {round(1:sk(1):df(1)), round(1:sk(2):df(2)), round(1:sk(3):df(3))};
%==========================================================================

%==========================================================================
function [ll,mf,vf] = inu_recon(f,chan,T,Sig)
% INU reconstruction
d  = [size(chan(1).B1,1) size(chan(1).B2,1) size(chan(1).B3,1)];
nz = d(3);
C  = numel(T);
if nargin<4,  Sig = cell(1,C); end
if nargout>1, mf  = zeros([d C],'single'); end
if nargout>2, vf  = zeros([d C],'single'); end
ll = zeros(2,C);

for c=1:C
    t  = T{c};
    dt = [size(t,1) size(t,2) size(t,3)];
    S  = Sig{c};
    if ~isempty(t)
        B1      = chan(c).B1;
        B2      = chan(c).B2;
        B3      = chan(c).B3;
        ICO     = chan(c).ICO;

        % ln p(t|IC0). Note that there's no regularisation for the DC component, so to avoid dealing with
        % a singular matrix, only ICO(2:end,2:end) is used for computing the determinant.
        ll(1,c) = -0.5*t(:)'*ICO*t(:) + sum(log(diag(chol(ICO(2:end,2:end)))),'double') - 0.5*(size(ICO,1)-1)*log(2*pi);

        %% Use BIC to account for having no regularisation on the DC component.
        %% Need to check that this is the right thing to do.
        %ll(1,c) = ll(1,c) - 0.5*log(sum(sum(sum(isfinite(f(:,:,:,c))))));

        if ~isempty(S)
            approx = 0;
            if size(S,2)==1   % S is the diagonal to save memory
                s      = S;
                approx = 0;
            end

            if approx<0       % Full computation
                C       = chol(S);
                % -E[ln q(t)]: 0.5*log(det(S)) + 0.5*D*log(2*pi) + 0.5*trace(S\S)
                ll(1,c) = ll(1,c) + sum(log(diag(C)),'double') + 0.5*size(C,1) + 0.5*size(C,1)*log(2*pi);
                % make it E[ln p(t|ICO)] instead of just p(t|ICO)
                ll(1,c) = ll(1,c) - 0.5*sum(sum(ICO.*S,'double'),'double'); 
                C       = reshape(C',[dt size(C,1)]);
            else
               if approx==0  % Diagonal approximation
                    % -E[ln q(t)]
                    if size(S,2) == size(S,1)
                        s   = diag(S);
                    end
                    ll(1,c) = ll(1,c) + 0.5*sum(log(s),'double') + 0.5*size(s,1) + 0.5*size(s,1)*log(2*pi);
                    ll(1,c) = ll(1,c) - 0.5*sum(diag(ICO).*s,'double');         % make it E[ln p(t|ICO)]
                    U       = zeros([0 0 0 0]);
                    s       = reshape(s,dt);
                else
                    % -E[ln q(t)]
                    ll(1,c) = ll(1,c) + sum(log(diag(chol(S))),'double') + 0.5*size(S,1) + 0.5*size(S,1)*log(2*pi);
                    ll(1,c) = ll(1,c) - 0.5*sum(sum(ICO.*S,'double'),'double'); % make it E[ln p(t|ICO)]

                    % Approximate with a few eigenmodes + diagonal
                    approx = min(approx,size(S,2)-1);
                    [V,D]  = eigs(double(S),approx);
                    U      = single(V*diag(sqrt(diag(D))));
                    s      = reshape(single(diag(S - U*U')),dt);
                    U      = reshape(U,[dt approx]); 
                end

                % Basis functions for diagonal (part of the ) approximation
                B1s = B1.^2;
                B2s = B2.^2;
                B3s = B3.^2;
            end

            for z=1:nz % Loop over slices
 
                %% Symbolic workings for E[f.*exp(x)] & Var[f.*exp(x)]
                %  syms m x em1 real
                %  syms f S positive
                %  ex = simplify(int(exp(x)*exp(-0.5*(x-m)^2/S)/sqrt(2*pi*S),x,-Inf,Inf),1000)
                %% E[f.*exp(x)] = f.*exp(m+S/2)
                %  v  = simplify(int(f^2*(exp(x)-ex)^2*exp(-0.5*(x-m)^2/S)/sqrt(2*pi*S),x,-Inf,Inf),1000)
                %% Var[f.*exp(x)] = f.^2.*exp(2*m).*exp(S).*(exp(S)-1)

                if approx<0
                    vl = single(0);
                    for ii=1:size(C,4)
                        vl = vl + inu_transform(B1,B2,B3(z,:),C(:,:,:,ii)).^2;
                    end
                else
                    % Approximate voxl-wise variance estimates from covariance matrix S 
                    vl = inu_transform(B1s,B2s,B3s(z,:),s); % Diagonal approximation
                    for ii=1:size(U,4)
                        vl = vl + inu_transform(B1,B2,B3(z,:),U(:,:,:,ii)).^2;
                    end
                end

                ml         = inu_transform(B1,B2,B3(z,:),t);
                fz         = f(:,:,z,c);
                ll(2,c)    = ll(2,c) + sum(ml(isfinite(fz(:))),'double');
                if nargout>1, mf(:,:,z,c) = f(:,:,z,c).*(exp(ml+vl/2)); end
                if nargout>2, vf(:,:,z,c) = f(:,:,z,c).^2.*exp(2*ml).*exp(vl).*(exp(vl)-1); end
            end
        else
            for z=1:nz % Loop over slices
                fz         = f(:,:,z,c);
                ml         = inu_transform(B1,B2,B3(z,:),t);
                ll(2,c)    = ll(2,c) + sum(ml(isfinite(fz(:))),'double');
                if nargout>1, mf(:,:,z,c) = f(:,:,z,c).*exp(ml); end
            end
        end
    else
        if nargout>1, mf(:,:,:,c) = f(:,:,:,c); end
    end
end
%==========================================================================

%==========================================================================
function Xo = insert2cell(Xo,C,X,c)
codes = unique(C);
codes = codes(codes ~= 0);
for i=1:numel(codes)
    msk    = (C == codes(i));
    io     = code2bin(codes(i),size(X,2));
    if any(io==c)
        Xo{i}(:,io==c) = X(msk);
    end
end
%==========================================================================

%==========================================================================
function bin = code2bin(code, length)
% FORMAT bin = code2bin(code, length)
%
% Convert a "missing code" to a mask of observed channels

base = uint64(2).^uint64(0:(length-1));
bin  = bitand(uint64(code),base) > 0;


