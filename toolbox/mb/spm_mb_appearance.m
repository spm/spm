function varargout = spm_mb_appearance(varargin)
%__________________________________________________________________________
%
% Functions for appearance model related.
%
% FORMAT chan       = spm_mb_appearance('BiasBasis',T,df,Mat,reg,samp)
% FORMAT [inu,ll]    = spm_mb_appearance('BiasField',T,chan,d,varargin)
% FORMAT labels     = spm_mb_appearance('GetLabels',dat,sett,do_samp)
% FORMAT [dat,sett] = spm_mb_appearance('IntroduceMG',dat,sett)
% FORMAT z          = spm_mb_appearance('Responsibility',m,b,W,n,f,mu,msk_chn)
% FORMAT [z,dat]    = spm_mb_appearance('Update',dat,mu0,sett)
% FORMAT dat        = spm_mb_appearance('UpdatePrior',dat,sett)
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% $Id$

if nargin == 0
    help spm_mb_appearance
    error('Not enough argument. Type ''help spm_mb_appearance'' for help.');
end
id       = varargin{1};
varargin = varargin(2:end);
switch id
    case 'BiasBasis'
        [varargout{1:nargout}] = BiasBasis(varargin{:});
    case 'BiasField'
        [varargout{1:nargout}] = BiasField(varargin{:});
    case 'GetLabels'
        [varargout{1:nargout}] = GetLabels(varargin{:});
    case 'IntroduceMG'
        [varargout{1:nargout}] = IntroduceMG(varargin{:});
    case 'Responsibility'
        [varargout{1:nargout}] = Responsibility(varargin{:});
    case 'Update'
        [varargout{1:nargout}] = Update(varargin{:});
    case 'UpdatePrior'
        [varargout{1:nargout}] = UpdatePrior(varargin{:});
    case 'UpdateSharedPrior'
        [varargout{1:nargout}] = UpdateSharedPrior(varargin{:});
    otherwise
        help spm_mb_appearance
        error('Unknown function %s. Type ''help spm_mb_appearance'' for help.', id)
end
end
%==========================================================================

%==========================================================================
% BiasField()
function [inu,ll] = BiasField(T,chan)
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
        inu_c        = TransformBF(chan(c).B1,chan(c).B2,chan(c).B3(z,:),T{c});
        inu(:,:,z,c) = single(exp(inu_c));
    end
end
end
%==========================================================================

%==========================================================================
% TransformBF()
function t = TransformBF(B1,B2,B3,T)
% Create an image-space log bias field from its basis function encoding.
%
% FORMAT t = TransformBF(B1,B2,B3,T)
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
end
%==========================================================================

%==========================================================================
function chan = BiasBasis(T,df,Mat,reg,samp)
if nargin<5, samp = 0; end
if nargin<4, reg  = ones(1,numel(T)); end
cl   = cell(1, numel(T));
chan = struct('ICO', cl, 'B1',cl, 'B2',cl, 'B3',cl);
ind  = SampleInd(df,samp);
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
end
%==========================================================================

%==========================================================================
% GetLabels()
function labels = GetLabels(dat,sett,samp)
if nargin < 3, samp = [1 1 1]; end
lab = dat.model.gmm.lab;
if isempty(lab), labels = 0; return; end

cm_map     = lab.cm_map; % cell array that defines the confusion matrix (cm)

% Load labels
labels = spm_mb_io('GetData',lab.f);
labels = SubSample(labels,samp);

% Use labels2use to keep only labels of interest
labels = round(labels(:));
labels(labels<1 || labels>numel(cm_map)) = numel(cm_map)+1;

% Get confusion matrix that maps from label value to probability value
cm     = GetLabelConfMatrix(cm_map,sett);

% Build NxK1 label image using confusion matrix
labels = cm(labels,:);

end
%==========================================================================

%==========================================================================
% GetLabelConfMatrix()
function cm = GetLabelConfMatrix(cm_map,w,K)
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
cm = zeros([L+1 K1],'single'); % Allocate confusion matrix
for l=1:L % Loop over labels
    ix            = false(1,K1);
    ix(cm_map{l}) = true;
    cm(l, ix)     = log(w/sum( ix));
    cm(l,~ix)     = log((1-w)/sum(~ix));
end
cm(L+1,:) = zeros(1,K1);
end
%==========================================================================

%==========================================================================
% Responsibility()
function [z,lb] = Responsibility(m,b,W,n,f,mu,msk_chn)
% Compute responsibilities.
%
% FORMAT z = Responsibility(m,b,W,n,f,mu,L,code)
% m       - GMM Means
% b       - GMM Mean d.f.
% W       - GMM Scale matrices
% n       - GMM Scale d.f.
% f       - Bias-corrected observed image in matrix form [nbvox nbchannel]
% mu      - Deformed and exponentiated template
% msk_chn - Mask of observed channels per code
% z       - Image of responsibilities [nbvox K]

const  = spm_gmm_lib('Normalisation', {m,b}, {W,n}, msk_chn);
f      = spm_gmm_lib('Marginal', f, {m,W,n}, const, msk_chn);
[z,lb] = spm_gmm_lib('Responsibility', f, mu);
end
%==========================================================================

%==========================================================================
function [dat,Z] = Update(dat,mu0,sett)
% Update appearance model for a single subject (GMM & bias field)
%
% FORMAT [dat,Z] = Update(dat,mu0,sett)
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
nit_lsinu    = 8;
reg          = gmm.inu_reg;
scal_inu     = 0.99;
K            = size(mu0,4);
Kmg          = numel(mg_ix);
Mat          = dat.Mat;
df           = dat.dm;
ds           = [size(mu0) 1 1];
ds           = ds(1:3);

% Get image data
f0    = spm_mb_io('GetImage',gmm);
f0    = SubSample(f0,dat.samp);
samp  = dat.samp;
samp2 = dat.samp2;

% GMM posterior: initial estimates at each iteration are made to be uncertain
m    = gmm.m;
b    = gmm.b*0+1e-3;
nval = C-1+1e-3;
%W   = repmat(mean(bsxfun(@times,gmm.W,reshape(gmm.n,[1 1 Kmg])),3)/nval,[1 1 Kmg]);
W    = bsxfun(@times,gmm.W,reshape(gmm.n,[1 1 numel(gmm.n)])./nval);
n    = gmm.n*0+nval;
mg_w = gmm.mg_w;


% Things to test
if true
    % If template is missing, set corresponding voxels of
    % f to missing too.
    msk_vx = any(isnan(mu0),4);
    for c=1:size(f0,4)
        fc = f0(:,:,:,c);
        fc(msk_vx) = NaN;
        f0(:,:,:,c) = fc;
    end
   %scl_samp = numel(msk_vx)/(numel(msk_vx)-sum(msk_vx(:)));
    clear msk_vx
else
    % Assume template is zero outside field of view
    mu0(~isfinite(mu0)) = 0;
   %scl_samp = 1.0;
end

% Compute Gaussian parameters on a subset of the voxels
[f,d] = SubSample( f0,samp2);
mu    = SubSample(mu0,samp2);

% Template
mu    = vol2vec(spm_mb_shape('TemplateK1',mu)); % Make K + 1 template
mu    = mu + GetLabels(dat,sett,true); % Add labels and template
mu    = mu(:,mg_ix); % Expand, if using multiple Gaussians per tissue

% Bias field related
do_inu = ~cellfun(@isempty,gmm.T);
if any(do_inu)
    T            = gmm.T;
    chan         = BiasBasis(T,df,Mat,reg,samp.*samp2);
    [inu,llpinu] = BiasField(T,chan);
    inuf         = inu.*f;
else
    inuf         = f;
end

% Format for spm_gmm
[inuf,code_image,msk_chn] = spm_gmm_lib('obs2cell', vol2vec(inuf));
mu                        = spm_gmm_lib('obs2cell', mu, code_image, false);
code_list                 = unique(code_image);
code_list                 = code_list(code_list ~= 0);

lbs = -Inf;
if any(do_inu)
    % Make sure bias field part of lower bound is correct
    lxb = SumLnINU(vol2vec(inu),code_image) + sum(llpinu);
else
    lxb = 0;
end

for it_appear=1:nit_appear

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update GMM and get responsibilities (Z)
    %------------------------------------------------------------
    lbso = lbs;
    pr   = sett.gmm(gmm.pop).pr;
    [Z,mog,~,lb,mg_w] = spm_gmm_lib('loop',inuf,1,...
                                    {{m,b},{W,n}},{'LogProp', mu}, ...
                                    'GaussPrior',   pr, ...
                                    'Missing',      msk_chn, ...
                                    'IterMax',      nit_gmm, ...
                                    'Tolerance',    tol_gmm, ...
                                    'SubIterMax',   nit_gmm_miss, ...
                                    'SubTolerance', tol_gmm, ...
                                    'Verbose',      [0 0], ...
                                    'MultGaussPi',  {mg_ix,mg_w});
    m      = mog.MU;
    b      = mog.b;
    W      = mog.V;
    n      = mog.n;
    lx     = lb.X(end)+lb.Z(end)+lb.P(end);
    lbs    = lb.sum(end)+lxb;
    if (it_appear==nit_appear) || (abs(lbs - lbso)/abs(lbs) < tol_gmm) 
        % Finished
        break
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update bias field parameters
    % This computes the derivatives of the negative logarithm of the
    % joint probability distribution
    %------------------------------------------------------------
    if any(do_inu)

        % Update bias field parameters for each channel separately
        for c=1:C % Loop over channels
            if isempty(T{c}), continue; end

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
                    Woo = W(ixo,ixo,k);
                    Wom = W(ixo,ixm,k);
                    Wmm = W(ixm,ixm,k);
                    Wmo = W(ixm,ixo,k);
                    Ao  = Woo - Wom*(Wmm\Wmo);
                    Ao  = (n(k) - nm) * Ao;
                    mo  = m(ixo,k);

                    % Compute statistics
                    gk = bsxfun(@minus, inuf{l}, mo.') * Ao(ixc,:).';
                    Hk = Ao(ixc,ixc);

                    gk = bsxfun(@times, gk, Z{l}(:,k));
                    Hk = bsxfun(@times, Hk, Z{l}(:,k));

                    % Accumulate across clusters
                    go = go + gk;
                    Ho = Ho + Hk;
                end

                % Multiply with bias corrected value (chain rule)
                oinuf = inuf{l}(:,ixc);
                go    = go .*  oinuf;
                Ho    = Ho .* (oinuf.^2);
                clear oinuf

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
            clear b3

            % Gauss-Newton update of bias field parameters
            Update = reshape((H + chan(c).ICO)\(gr + chan(c).ICO*T{c}(:)),size(T{c}));
            clear H gr

            % Line-search
            armijo  = scal_inu;
            oT      = T{c};
            ollpinu = llpinu;
            olx     = lx;
            olxb    = lxb; 
            for ls=1:nit_lsinu

                % Update bias-field parameters
                T{c} = T{c} - armijo*Update;

                % Compute new bias-field (only for channel c)
                [inu(:,:,:,c),llpinu(c)] = BiasField(T(c),chan(c));
                inuf  = spm_gmm_lib('obs2cell', vol2vec(inu.*f), code_image, true);

                % Recompute responsibilities and lower bounds
                % (with updated bias field)
                [Z,lx] = Responsibility(m,b,W,n,inuf,ReWeightMu(mu,log(mg_w)),msk_chn);
                lxb    = SumLnINU(vol2vec(inu), code_image) + sum(llpinu);
                % Check new lower bound
                if  ((lx + lxb) - (olx + olxb))/abs(lx + lxb) > -eps('single')
                    % Converged
                    break;
                else
                    armijo = armijo*0.5;
                    T{c}   = oT;
                    if ls == nit_lsinu
                        % Did not converge -> reset
                        llpinu       = ollpinu;
                        inu(:,:,:,c) = BiasField(T(c),chan(c));
                        inuf         = spm_gmm_lib('obs2cell', vol2vec(inu.*f), code_image, true);

                        % Recompute responsibilities
                        Z = Responsibility(m,b,W,n,inuf,ReWeightMu(mu,log(mg_w)),msk_chn);
                    end
                end
            end
            clear oT Update
        end
        gmm.T = T;
    end
end
clear f inu mu

% Update dat
lbs      = lx+lxb+lb.MU(end)+lb.A(end);
gmm.m    = m;
gmm.b    = b;
gmm.W    = W;
gmm.n    = n;
gmm.lb   = lb;
gmm.mg_w = mg_w;
dat.model.gmm = gmm;

if nargout > 1
    % Compute full size resps

    if any(samp2~=1)
        % Compute full-sized responsibilities on original data

        if any(do_inu == true) % Bias correct
            chan = BiasBasis(T,df,Mat,reg,samp);
            [inu,llpinu] = BiasField(T,chan);
            lxb  = SumLnINU(vol2vec(inu), code_image) + sum(llpinu);
            inuf = inu.*f0;
        else
            inuf = f0;
        end

        [inuf,code_image,msk_chn] = spm_gmm_lib('obs2cell', vol2vec(inuf));

        % Template
        mu0 = vol2vec(spm_mb_shape('TemplateK1',mu0));
        mu0 = mu0 + GetLabels(dat,sett);
       %mu0 = mu0(:,mg_ix) + log(mg_w);
        mu0 = bsxfun(@plus, mu0(:,mg_ix), log(mg_w));
        mu0 = spm_gmm_lib('obs2cell', mu0, code_image, false);

        [Z,lx]   = Responsibility(m,b,W,n,inuf,mu0,msk_chn);
        lbs      = lx+lxb+lb.MU(end)+lb.A(end);
        clear mu0
    end
    Z = spm_gmm_lib('cell2obs', Z, code_image, msk_chn);

    % If using multiple Gaussians per tissue, collapse so that Z is of
    % size K.
    % Note: This assumes that mg_ix is sorted into increasing values.
    if size(Z,2) > K
        for k=1:K
            Z(:,k) = sum(Z(:,mg_ix==k),2);
        end
    end
    Z = Z(:,1:K);

    % Make 4D
    Z = vec2vol(Z,ds);
end
dat.E(1) = -lbs; %*scl_samp;
end
%==========================================================================

%==========================================================================
function X2d = vol2vec(X4d)
d   = [size(X4d) 1 1];
X2d = reshape(X4d,[prod(d(1:3)) d(4)]);
end
%==========================================================================

%==========================================================================
function X4d = vec2vol(X2d,dm)
dm = [dm(:)' 1 1];
if size(X2d,1)~=prod(dm(1:3))
    error('Incompatible dimensions.');
end
X4d = reshape(X2d,[dm(1:3) size(X2d,2)]);
end
%==========================================================================

%==========================================================================
function sett = UpdatePrior(dat, sett)
if isempty(sett.gmm), return; end

% Get population indices
code = zeros(numel(dat),1);
for n=1:numel(dat)
    if isfield(dat(n).model,'gmm')
        code(n) = dat(n).model.gmm.pop;
    end
end

for p=1:numel(sett.gmm) % Loop over populations
    if sett.gmm(p).updt_intens
        index = find(code==p);
        N     = numel(index);
        pr    = sett.gmm(p).pr;

        % Get all posteriors
        po    = cell(1,N);
        for n=1:N
            n1          = index(n);
            po{n}{1}{1} = dat(n1).model.gmm.m;
            po{n}{1}{2} = dat(n1).model.gmm.b;
            po{n}{2}{1} = dat(n1).model.gmm.W;
            po{n}{2}{2} = dat(n1).model.gmm.n;
        end

        % Update prior
        sett.gmm(p).pr = spm_gmm_lib('updatehyperpars',po,pr);
    end
end
end
%==========================================================================

%==========================================================================
function sett = UpdateSharedPrior(dat, sett)
if isempty(sett.gmm), return; end

% Get population indices
code = zeros(numel(dat),1);
for n=1:numel(dat)
    if isfield(dat(n).model,'gmm')
        code(n) = dat(n).model.gmm.pop;
    end
end

for p=1:numel(sett.gmm) % Loop over populations
    if sett.gmm(p).updt_intens
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
                po{nk} = {{gmm.m(:,k),gmm.b(:,k)},{diag(diag(gmm.W(:,:,k))),gmm.n(:,k)}};
            end
        end

        % Update prior
        pr    = spm_gmm_lib('updatehyperpars',po,pr);
        pr{1} = repmat(pr{1},1,K);
        pr{2} = repmat(pr{2},1,K);
        pr{3} = repmat(pr{3},1,1,K);
        pr{4} = repmat(pr{4},1,K);
        sett.gmm(p).pr = pr;
    end
end
end
%==========================================================================

%==========================================================================
%
% Utility functions
%
%==========================================================================

%==========================================================================
function lb = SumLnINU(varargin)
% Compute parts of the lower bound
%
% FORMAT lb = SumLnINU(inu,code_image)
% inu        - Exponentiated INU field [one channel]
% code_image - integer type data indicating whether data present
% lb         - Sum of the log bias field
%   >> This is part of the normalisation term in the GMM
%
inu        = varargin{1};
code_image = varargin{2};
lb         = 0;
if islogical(code_image)
    for c=1:size(inu,2)
        lb   = lb + sum(log(inu(code_image,c)),'double');
    end
else
    for c=1:size(inu,2)
        code = zeros(1,'like',code_image)+2^(c-1);
        msk  = bsxfun(@bitand,code_image,code)>0;
        lb   = lb + sum(log(inu(msk,c)),'double');
    end
end
end
%==========================================================================

%==========================================================================
function mu = ReWeightMu(mu,logmg_w)
if sum(logmg_w) == 0, return; end
for i=1:numel(mu)
   %mu{i} = mu{i} + logmg_w;
    mu{i} = bsxfun(@plus, mu{i}, logmg_w);
end
end
%==========================================================================

%==========================================================================
function [of,d,w] = SubSample(f,samp)
% Subsample a multichannel volume.
%
% FORMAT [of,d,scl_samp] = SubSample(f,samp);
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

    ind  = SampleInd(df,samp);
    d    = cellfun(@length,ind);  % New dimensions
    of   = f(ind{:},:,:);

    % For weighting data parts of lowerbound with factor based on amount of
    % downsampling
    w    = prod(df)/prod(d);
end
end
%==========================================================================

%==========================================================================
function ind = SampleInd(df,samp)
sk   = max([1 1 1],samp);
ind  = {round(1:sk(1):df(1)), round(1:sk(2):df(2)), round(1:sk(3):df(3))};
end
%==========================================================================

