function varargout = spm_mb_appearance(varargin) % Appearance model
%
% FORMAT chan       = spm_mb_appearance('inu_basis',T,df,Mat,reg,samp)
% FORMAT [inu,ll]   = spm_mb_appearance('inu_field',T,chan,d,varargin)
% FORMAT z          = spm_mb_appearance('responsibility',m,b,W,n,f,mu,msk_chn)
% FORMAT dat        = spm_mb_appearance('restart',dat,sett)
% FORMAT [z,dat]    = spm_mb_appearance('update',dat,mu,sett)
% FORMAT dat        = spm_mb_appearance('update_prior',dat,sett)
% FORMAT            = spm_mb_appearance('debug_show',img,img_is,modality,fig_title,do)
%__________________________________________________________________________

% Mikael Brudfors, John Ashburner & Yael Balbastre
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


[varargout{1:nargout}] = spm_subfun(localfunctions,varargin{:});
%==========================================================================

%==========================================================================
function [dat,sett] = restart(dat,sett)
% Allow a restart of the GMM
for n=1:numel(dat)
    if isfield(dat(n).model,'gmm')
        p    = dat(n).model.gmm.pop;
        same = ~all(sum(diff(sett.gmm(p).pr{1},1,2).^2,1));
        if same
            % Intensity priors are identical for all clusters
            % so need to break the symmetry.
            % Increase variance of Gaussians and make the means
            % more similar to each other.
            dat(n).model.gmm.W = dat(n).model.gmm.W*0.1;
            dat(n).model.gmm.m = bsxfun(@plus, dat(n).model.gmm.m*0.1,...
                                        mean(dat(n).model.gmm.m,2)*0.9);
            dat(n).model.gmm.Alpha = dat(n).model.gmm.Alpha*0+eps;
        else
            if false % numel(sett.gmm)==1
                % Initialise to the current priors.
                % Tissue priors have less effect on the posterior
                % whereas intensity priors have a stronger effect
                % on ensuring a consistent clustering across subjects.
                % Commented out because it appears to cause strange
                % behaviour for empty clusters.
                dat(n).model.gmm.m  = sett.gmm(p).pr{1};
                dat(n).model.gmm.b  = sett.gmm(p).pr{2};
                dat(n).model.gmm.W  = sett.gmm(p).pr{3};
                dat(n).model.gmm.nu = sett.gmm(p).pr{4};
                dat(n).model.gmm.Alpha = dat(n).model.gmm.Alpha*0+eps;
           else
                % Set all GMM parameters to be identical so the
                % next GMM fit restarts everything at the current
                % priors. Symmetry already broken with intensity
                % priors.
                % Tissue priors have more effect on the posterior
                % which reduces the chance of different populations
                % having their own clustering. This is an attempt
                % to ensure that clusters correspond across
                % populations.
                [M,K1] = size(dat(n).model.gmm.m);
                dat(n).model.gmm.m  = ones(M,K1)*1000;
                dat(n).model.gmm.b  = ones(1,K1);
                dat(n).model.gmm.W  = repmat(eye(M,M)/1000,[1 1 K1]);
                dat(n).model.gmm.nu = ones(1,K1)*(M+1);
                dat(n).model.gmm.Alpha = dat(n).model.gmm.Alpha*0+eps;
           end
        end
    end
end
%==========================================================================

%==========================================================================
function [Z,lx] = responsibility(m,b,W,nu, mg_ix, mu,mf,vf)
% Compute responsibilities.
%
% FORMAT z = responsibility(m,b,W,nu, mg_ix, mu,mf,vf)
% m       - GMM Means
% b       - GMM Mean d.f.
% W       - GMM Scale matrices
% n       - GMM Scale d.f.
% mg_ix   - GMM indices
% mu      - Deformed and template
% mf      - Expectation of INU-corrected observed image in matrix form [nbvox nbchannel]
% vf      - Variance of INU-corrected observed image in matrix form [nbvox nbchannel]
% Z       - Image of responsibilities [nbvox K]
[Z,lx] = spm_gmmlib('resp',m,b,W,nu, uint64(mg_ix),mu, mf,vf, uint64([1 1 1]));
%==========================================================================

%==========================================================================
function [dat,Z] = update(dat,mu,sett)
% Update appearance model for a single subject (GMM & bias field)
%
% FORMAT [dat,Z] = update(dat,mu,sett)
% dat - Structure holding data for a single subject
% mu  - Log template
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
Mat          = dat.Mat;
df           = dat.dm;
ds           = [size(mu) 1 1];
ds           = ds(1:3);

% For visual debugging (disable/enable in debug_show())
debug_show(mu, 'template_k1');

% Get image data
samp1 = dat.samp;
samp2 = gmm.samp;
f0    = subsample(spm_mb_io('get_image',gmm),samp1);

if isa(dat.lab,'struct')
    label = spm_mb_io('get_data', dat.lab.f);
    label = uint8(subsample(label,samp1));
else
    label = [];
end

% For visual debugging (disable/enable in debug_show())
debug_show(f0, 'observed', gmm.modality);

% Intensity priors
pr      = sett.gmm(gmm.pop).pr;

% GMM posterior
cluster = {gmm.m,gmm.b,gmm.W,gmm.nu,gmm.gam,gmm.Alpha};

% Bias field related
T     = gmm.T;
if isfield(gmm,'Sig')
    Sig = gmm.Sig;
else
    Sig = cell(size(T));
end

do_inu     = ~cellfun(@isempty,T);
ind        = sample_ind(ds,samp2);
[msk,nvox] = data_mask(f0,mu,samp2);
if any(do_inu)
    chan          = inu_basis(T,df,Mat,inu_reg,samp1.*samp2);
    [llinu,mf,vf] = inu_recon(f0(ind{:},:),msk,chan,T,Sig);
    lxb           = sum(llinu(:),'double');
else
    mf            = f0(ind{:},:);
    vf            = zeros(size(mf),'single');
    lxb           = 0;
end
lbs  = -Inf;

for it_appear=1:nit_appear
    [cluster,lb] = spm_mb_gmm(mf,vf, cluster, pr, uint64(mg_ix), mu, label, ...
                             samp2, nit_gmm, nvox*tol_gmm, nit_gmm_miss, nvox*tol_gmm*0.1);
    lnP    = bsxfun(@minus, psi(cluster{6}), psi(sum(cluster{6},1)));
    if ~isempty(lnP), lnP(1,:) = 0; end

    lx     = lb.X+lb.P+lb.Alpha;
    lbso   = lbs;
    lbs    = lb.mu+lb.A+lx+lxb;
   %fprintf(' %8.4f', -lbs/nvox);

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
        for c=1:C % Loop over channels

            if isempty(T{c}), continue; end

            % Compute gradient and Hessian
            d3 = numel(T{c}); % Number of DCT parameters
            H  = zeros(d3,d3);
            gr = zeros(d3,1);

            B1 = double(chan(c).B1);
            B2 = double(chan(c).B2);
            B3 = double(chan(c).B3);

            lx   = 0;
            for z=1:size(mf,3)
                if isempty(label)
                    labz = [];
                else
                    labz = label(ind{1:2},ind{3}(z));
                end
                [gr_im,H_im,lxt] = spm_gmmlib('inugrads', cluster{1:5}, uint64(mg_ix),mu(ind{1:2},ind{3}(z),:),...
                                              mf(:,:,z,:),vf(:,:,z,:), uint64([1 1 1]), labz,lnP, uint64(c));
                lx = lx + lxt;
                b3 = B3(z,:)';
                gr = gr + kron(b3,spm_krutil(double(gr_im),B1,B2,0));
                H  = H  + kron(b3*b3',spm_krutil(double(H_im),B1,B2,1));
            end
            clear gr_im H_im b3

            % Gauss-Newton update of bias field parameters
            Sig{c} = 1./diag(H + chan(c).L); % Diagonal approximation
            T{c}   = T{c} - reshape((H + chan(c).L)\(gr + chan(c).L*T{c}(:)),size(T{c}));
            clear H gr

            % Compute new expectations (only for channel c)
            [llinu(:,c),mf(:,:,:,c),vf(:,:,:,c)] = inu_recon(f0(ind{:},c),msk,chan(c),T(c),Sig(c));
            lxb    = sum(llinu(:),'double');
        end
    end
end

% Update dat
lbs       = lx+lxb+lb.mu+lb.A+lb.Alpha;
gmm.T     = T;
gmm.Sig   = Sig;
gmm.m     = cluster{1};
gmm.b     = cluster{2};
gmm.W     = cluster{3};
gmm.nu    = cluster{4};
gmm.gam   = cluster{5};
gmm.Alpha = cluster{6};
gmm.lb    = lb;
if nargout > 1
    [msk,nvox] = data_mask(f0,mu);

    % Compute full-sized responsibilities on original data
    if any(do_inu) % Bias correct
        chan          = inu_basis(T,df,Mat,inu_reg,samp1);
        [llinu,mf,vf] = inu_recon(f0,msk,chan,T,Sig);
        lxb = sum(llinu(:));
    else
        mf  = f0;
        vf  = zeros(size(mf),'single');
        lxb = 0;
    end
    [Z,lx] = spm_gmmlib('resp',cluster{1:5}, uint64(mg_ix),mu, mf,vf, uint64(samp2), label,lnP);
    lbs    = lx+lxb+lb.mu+lb.A+lb.Alpha;
end

% For visual debugging (disable/enable in debug_show())
debug_show(Z,'responsibilities');

% Assign
dat.E(1)      = -lbs;
dat.nvox      = nvox;
dat.model.gmm = gmm;
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
        po    = cell(1,N);
        nlab  = 0;
        for n=1:N
            n1       = index(n);
            po{n}{1} = dat(n1).model.gmm.m;
            po{n}{2} = dat(n1).model.gmm.b;
            po{n}{3} = dat(n1).model.gmm.W;
            po{n}{4} = dat(n1).model.gmm.nu;
            if ~isempty(dat(n1).model.gmm.Alpha), nlab = nlab + 1; end
        end

        % Update prior
        hp = sett.gmm(p).hyperpriors;
        sett.gmm(p).pr = spm_mb_gmm('updatehyperpars',po,pr,hp{:});

        % Attempt to increase stability by avoiding singular precision matrices
        W  = sett.gmm(p).pr{3};
        for k=1:size(W,3)
            S        = inv(W(:,:,k));
            W(:,:,k) = inv(S*(1-1e-9) + 1e-9*mean(diag(S))*eye(size(S)));
        end
        sett.gmm(p).pr{3} = W;

        %% Update INU regularisation. Disabled because it under-regularises
        %ss_inu0 = zeros(1,size(pr{1},1));
        %ss_inu1 = zeros(1,size(pr{1},1));
        %for n=1:N
        %    [tr,D] = inu_reg_suffstat(dat(n1).model.gmm.T,dat(n1).model.gmm.Sigma,dat(n1).df,dat(n1).Mat);
        %    ss_inu1 = ss_inu1 + tr;
        %    ss_inu0 = ss_inu0 + D;
        %end
        %sett.gmm(p).inu_reg = ss_inu1./ss_inu2;

        % Dirichlet hyperpriors
        if nlab>0
            lb     = 0;
            K1     = size(sett.gmm(p).pr{1},2);
            Alpha  = zeros(255,nlab);
            Alpha0 = zeros(256,K1);
            for k=1:K1
                n2 = 1;
                for n=1:N
                    n1 = index(n);
                    if ~isempty(dat(n1).model.gmm.Alpha)
                        Alpha(:,n2) = dat(n1).model.gmm.Alpha(2:end,k);
                        n2          = n2 + 1;
                    end
                end
                [Alpha0(2:end,k),lbk] = spm_mb_gmm('dirichlet_hyperparameters',Alpha);
                lb = lb + lbk;
            end
            sett.gmm(p).pr{6} = Alpha0;
        else
            sett.gmm(p).pr{6} = [];
        end
    end
end
%==========================================================================

%==========================================================================
function [tr,D] = inu_reg_suffstat(T,Sig,df,Mat)
vs   = sqrt(sum(Mat(1:3,1:3).^2,1));
C    = numel(T);
tr   = zeros(1,C);
D    = zeros(1,C);
for c=1:numel(T)
    d3   = size(T{c});
    D(c) = prod(d3); % Should really subtract 1 because of the DC.
    if prod(d3)>0
        L = inu_regulariser(df,d3,vs);
        t = T{c}(:);
        S = Sig{c};
        if size(S,2)==1
            % Diagonal approximation
            tr(c) = sum(S.*t.^2) + sum(S.*diag(L));
        else
            tr(c) = t'*L*t + sum(sum(S.*L));
        end
    end
end
%==========================================================================

%==========================================================================
function [msk,nvox] = data_mask(mf,mu,samp2)
if nargin>=3
    ind = sample_ind(size(mf),samp2);
else
    ind = {':',':',':'};
end
msk  = isfinite(mu(ind{:},1));
if nargout>=2
    mskf = false(size(msk));
    for i=1:size(mf,4)
        mskf = mskf | isfinite(mf(ind{:},i));
    end
    nvox = sum(mskf(:) & msk(:));
end
%==========================================================================

%==========================================================================
function [of,d] = subsample(f,samp)
% Subsample a multichannel volume.
%
% FORMAT [of,d,scl_samp] = subsample(f,samp);
% f    - Original volume
% samp - Sampling distances in voxels
% of   - Resampled volume
% d    - Output dimensions

if all(samp==1)
    of = f;
    d  = [size(f,1) size(f,2) size(f,3)];
else
    % Input image properties
    df   = [size(f) 1];
    df   = df(1:3);

    ind  = sample_ind(df,samp);
    d    = cellfun(@length,ind);  % New dimensions
    of   = f(ind{:},:,:);
end
%==========================================================================

%==========================================================================
function ind = sample_ind(df,samp)
% Indices for sampling
% FORMAT ind = sample_ind(df,samp)
% df   - image dimensions
% samp - spacing along the three dimensions
% ind  - indices
df   = [df(:)' 1 1];
df   = df(1:3);
samp = [samp(:)' 1 1];
samp = samp(1:3);
sk   = max([1 1 1],samp);
ind  = {round(1:sk(1):df(1)), round(1:sk(2):df(2)), round(1:sk(3):df(3))};
%==========================================================================

%==========================================================================
function chan = inu_basis(T,df,Mat,reg,samp)
if nargin<5, samp = 0; end
if nargin<4, reg  = ones(1,numel(T)); end
cl   = cell(1, numel(T));
chan = struct('L', cl, 'B1',cl, 'B2',cl, 'B3',cl);
ind  = sample_ind(df,samp);
vs   = sqrt(sum(Mat(1:3,1:3).^2,1));
for c=1:numel(T)
    d3 = [size(T{c}) 1];
    d3 = d3(1:3);

    % Regulariser
    chan(c).L  = reg(c)*inu_regulariser(df,d3,vs);

    % Basis functions for bias correction
    chan(c).B1 = spm_dctmtx(df(1),d3(1),ind{1});
    chan(c).B2 = spm_dctmtx(df(2),d3(2),ind{2});
    chan(c).B3 = spm_dctmtx(df(3),d3(3),ind{3});
end
%==========================================================================

%==========================================================================
function L = inu_regulariser(df,d3,vs)
% Bending energy regulariser (precision matrix)
kx = (pi*((1:d3(1))'-1)/df(1)/vs(1)).^2;
ky = (pi*((1:d3(2))'-1)/df(2)/vs(2)).^2;
kz = (pi*((1:d3(3))'-1)/df(3)/vs(3)).^2;
L  =  (kron(kz.^2,kron(ky.^0,kx.^0)) + kron(kz.^0,kron(ky.^2,kx.^0)) + kron(kz.^0,kron(ky.^0,kx.^2)) +...
    2*(kron(kz.^1,kron(ky.^1,kx.^0)) + kron(kz.^1,kron(ky.^0,kx.^1)) + kron(kz.^0,kron(ky.^1,kx.^1))) );
L  = diag(max(L,0));
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
    ll(c) = -0.5*t(:)'*chan(c).L*t(:);
    for z=1:nz
        inu_c        = inu_transform(chan(c).B1,chan(c).B2,chan(c).B3(z,:),T{c});
        inu(:,:,z,c) = single(exp(inu_c));
    end
end
%==========================================================================

%==========================================================================
function [ll,mf,vf] = inu_recon(f,msk,chan,T,Sig)
% INU reconstruction
d  = [size(chan(1).B1,1) size(chan(1).B2,1) size(chan(1).B3,1)];
nz = d(3);
C  = numel(T);
if nargin<5 || isempty(Sig), Sig = cell(1,C); end
if nargout>1, mf  = zeros([d C],'single'); end
if nargout>2, vf  = zeros([d C],'single'); end
ll = zeros(2,C);

for c=1:C
    t  = T{c};
    dt = [size(t,1) size(t,2) size(t,3)];
    S  = Sig{c};
    if ~isempty(t)
        B1 = chan(c).B1;
        B2 = chan(c).B2;
        B3 = chan(c).B3;
        L  = chan(c).L;

        % ln p(t|L). Note that there's no regularisation for the DC component, so to avoid dealing with
        % a singular matrix, only L(2:end,2:end) is used for computing the determinant.
        ll(1,c) = -0.5*t(:)'*L*t(:) + sum(log(diag(chol(L(2:end,2:end)))),'double') - 0.5*(size(L,1)-1)*log(2*pi);

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
                % make it E[ln p(t|L)] instead of just p(t|L)
                ll(1,c) = ll(1,c) - 0.5*sum(sum(L.*S,'double'),'double');
                C       = reshape(C',[dt size(C,1)]);
            else
               if approx==0  % Diagonal approximation
                    % -E[ln q(t)]
                    if size(S,2) == size(S,1)
                        s   = diag(S);
                    end
                    ll(1,c) = ll(1,c) + 0.5*sum(log(s),'double') + 0.5*size(s,1) + 0.5*size(s,1)*log(2*pi);
                    ll(1,c) = ll(1,c) - 0.5*sum(diag(L).*s,'double');         % make it E[ln p(t|L)]
                    U       = zeros([0 0 0 0]);
                    s       = reshape(s,dt);
                else
                    % -E[ln q(t)]
                    ll(1,c) = ll(1,c) + sum(log(diag(chol(S))),'double') + 0.5*size(S,1) + 0.5*size(S,1)*log(2*pi);
                    ll(1,c) = ll(1,c) - 0.5*sum(sum(L.*S,'double'),'double'); % make it E[ln p(t|L)]

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
                    % Approximate voxel-wise variance estimates from covariance matrix S
                    vl = inu_transform(B1s,B2s,B3s(z,:),s); % Diagonal approximation
                    for ii=1:size(U,4)
                        vl = vl + inu_transform(B1,B2,B3(z,:),U(:,:,:,ii)).^2;
                    end
                end

                ml         = inu_transform(B1,B2,B3(z,:),t);
                if ~isempty(msk)
                    mskz   = msk(:,:,z);
                else
                    mskz   = true(size(f,1),size(f,2));
                end
                fz         = f(:,:,z,c);
                ll(2,c)    = ll(2,c) + sum(ml(isfinite(fz(:)) & mskz(:)),'double');
                if nargout>1, mf(:,:,z,c) = fz.*(exp(ml + vl/2)); end
                if nargout>2, vf(:,:,z,c) = fz.^2.*exp(2*ml + vl).*(exp(vl) - 1); end
            end
        else
            for z=1:nz % Loop over slices
                if ~isempty(msk)
                    mskz   = msk(:,:,z);
                else
                    mskz   = true(size(f,1),size(f,2));
                end
                fz         = f(:,:,z,c);
                ml         = inu_transform(B1,B2,B3(z,:),t);
                ll(2,c)    = ll(2,c) + sum(ml(isfinite(fz(:)) & mskz(:)),'double');
                if nargout>1, mf(:,:,z,c) = fz.*exp(ml); end
            end
        end
    else
        if nargout>1, mf(:,:,:,c) = f(:,:,:,c); end
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
function debug_show(img,img_is,modality,fig_title,do_show)
% FORMAT debug_show(img,img_is,modality,fig_title,do_show)
%
% Show 4D image, can be (img_is):
% * 'observed' : Observed image data (if multi-channel, shows only first).
% * 'responsibilities' : Tissue responsibilities.
% * 'template_k1' : Tissue template with K + 1 classes.
% * 'template' : Tissue template with K classes.

if nargin < 2, img_is    = 'observed'; end
if nargin < 3, modality  = 1; end
if nargin < 4, fig_title = ''; end
if nargin < 5, do_show   = false; end
if ~do_show || ~any(strcmpi({'observed','responsibilities','template_k1','template'},img_is))
    return;
end
% Create/find figure
f = findobj('Type', 'Figure', 'Name', img_is);
if isempty(f)
    f = figure('Name', img_is, 'NumberTitle', 'off');
end
set(0, 'CurrentFigure', f);
clf(f);
% Image type specific
if strcmp(img_is,'responsibilities')
    img = cat(4,img,1 - sum(img,4));
elseif strcmp(img_is,'observed')
    c   = 1;  % Channel to show
    img = img(:,:,:,c);
elseif strcmp(img_is,'template')
    img = spm_mb_classes('template_k1',img,[]);
end
clim = [-Inf Inf];
if modality == 2
    clim = [1000 1100];  % CT scan
end
% What slice index (ix) to show
dm         = [size(img) 1 1 1];
ix         = round(0.5*dm(1:3));
% Colormap
if strcmp(img_is,'observed')
    colormap('gray')
else
    num_colors = dm(4) + 1;
    colormap(hsv(num_colors))
end
% Axis z
img1 = img(:,:,ix(3),:);
if ~strcmp(img_is,'observed')
    msk       = any(~isfinite(img1),4);
    [~,img1]  = max(img1,[],4);
    img1(msk) = num_colors;
end
subplot(131)
imagesc(img1,clim); axis off;
% Axis y
img2 = img(:,ix(2),:,:);
if ~strcmp(img_is,'observed')
    msk       = any(~isfinite(img2),4);
    [~,img2]  = max(img2,[],4);
    img2(msk) = num_colors;
end
img2 = squeeze(img2);
subplot(132)
imagesc(img2,clim); axis off;
title(fig_title)
% Axis x
img3 = img(ix(1),:,:,:);
if ~strcmp(img_is,'observed')
    msk       = any(~isfinite(img3),4);
    [~,img3]  = max(img3,[],4);
    img3(msk) = num_colors;
end
img3 = squeeze(img3);
subplot(133)
imagesc(img3,clim); axis off;
% Draw
drawnow
%==========================================================================
