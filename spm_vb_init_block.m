function [block] = spm_vb_init_block(Y,block)
% Initialise Variational Bayes for GLM-AR models
% FORMAT [block] = spm_vb_init_block(Y,block)
%
% Y      - [T x N] time series with T time points, N voxels
% block  - data structure (see spm_vb_glmar)
%__________________________________________________________________________
% Copyright (C) 2005-2014 Wellcome Trust Centre for Neuroimaging

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_init_block.m 6079 2014-06-30 18:25:37Z spm $

k = block.k;
p = block.p;
N = block.N;
T = block.T;
X = block.X;

% Default optimisation parameters
if ~isfield(block,'tol')
    block.tol = 0.0001;
end
if ~isfield(block,'maxits')
    block.maxits = 4;
end
if ~isfield(block,'verbose')
    block.verbose = 1;
end

if block.verbose
    disp('Initialising block');
    disp(' ');
end


%-Initialise approximate alpha posterior
%--------------------------------------------------------------------------
block.b_alpha = block.b_alpha_prior;
if block.update_alpha
    block.c_alpha = N/2 + block.c_alpha_prior;
else
    block.c_alpha = block.c_alpha_prior;
end
block.mean_alpha  = block.b_alpha.*block.c_alpha;

%-Initialise approximate beta posterior
%--------------------------------------------------------------------------
if strcmp(block.priors.A,'Discrete')
    S = block.priors.S;
    block.c_beta_prior = block.c_beta_prior*ones(1,S);
    block.b_beta_prior = block.b_beta_prior*ones(1,S);
    block.as = rand(block.p,S);

    block.b_beta     = block.b_beta_prior;
    if block.update_beta
        block.c_beta = 0.5*ones(block.p,1)*block.priors.N+block.c_beta_prior;
    else
        block.c_beta = block.c_beta_prior;
    end
else
    block.b_beta     = block.b_beta_prior;
    if block.update_beta
        block.c_beta = N/2 + block.c_beta_prior;
    else
        block.c_beta = block.c_beta_prior;
    end
end
block.mean_beta = block.b_beta.*block.c_beta;

%-Initialise approximate lambda posterior
%--------------------------------------------------------------------------
block.b_lambda = block.b_lambda_prior;
if block.update_lambda
    block.c_lambda = (T-block.p)/2 + block.c_lambda_prior;
else
    block.c_lambda = block.c_lambda_prior;
end
block.mean_lambda  = block.b_lambda.*block.c_lambda;

% Initialise approximate w posterior
%--------------------------------------------------------------------------
if block.verbose
    disp('Initialising regression coefficient posterior');
end
try
    Xp = block.Xp;
    X2 = block.X2;
catch
    [ux,dx,vx] = svd(X);
    ddx     = diag(dx);
    svd_tol = max(ddx)*eps*k;
    rank_X  = sum(ddx > svd_tol);
    ddxm    = diag(ones(rank_X,1)./ddx(1:rank_X));
    ddxm2   = diag(ones(rank_X,1)./(ddx(1:rank_X).^2));
    Xp      = vx(:,1:rank_X)*ddxm*ux(:,1:rank_X)';
    X2      = vx(:,1:rank_X)*ddxm2*vx(:,1:rank_X)';
end  

w_ols  = Xp*Y;
w_mean = w_ols;
Y_pred = X*w_mean;
v      = mean((Y-Y_pred).^2);
for n=1:N
    w_cov_temp     = v(n)*X2;
    block.w_cov{n} = w_cov_temp;
end
block.w_mean  = w_mean(:);
block.w_ols   = block.w_mean;
block.wk_mean = reshape(block.w_mean,k,N);
block.wk_ols  = reshape(block.w_ols,k,N);
% Initialize WGL
if isfield(block.Dw,'vxyz')
    block.Dw = spm_vb_spatial_precision ('Spatial - WGL',block.Dw.vxyz,block.wk_ols);
end
% Initialise AR coefficient posterior
if block.verbose
    disp('Initialising AR coefficients');
end
% Embed data
for pp=1:p
    dy(pp,1:T-p,:) = Y(p-pp+1:T-pp,:);
end
if p>0
    e = Y(p+1:T,:) - Y_pred(p+1:T,:);
    for n=1:N
        for pp=1:p
            if block.k > 1
                dyhat(pp,:) = w_mean(:,n)'*squeeze(block.dX(pp,:,:));
            else
                dyhat(pp,:) = w_mean(:,n)*squeeze(block.dX(pp,:,:))';
            end
        end
        E_tilde     = dy(:,:,n)-dyhat;
        iterm       = inv(E_tilde * E_tilde');
        block.ap_mean(:,n) = (iterm * E_tilde*e(:,n));  
        e_pred      = E_tilde' * block.ap_mean(:,n);
        v2          = mean((e(:,n) - e_pred).^2);
        block.a_cov{n} = v2 * iterm;
        block.a2{n} = block.ap_mean(:,n)*block.ap_mean(:,n)'+block.a_cov{n};
    end
    block.ap_ols    = block.ap_mean;
    block.a_mean    = block.ap_mean(:);
    
    if strcmp(block.priors.A,'block-Limiting')
        if p==1
            a_mean=mean(block.ap_mean);
            block.ap_mean=a_mean*ones(1,N);
        end
    end
end

%-Setting up cross-covariance matrices
%--------------------------------------------------------------------------
if p>0
    if block.verbose
        disp('Setting up cross-covariance matrices');
    end
    % Get input-output lagged covariances (I.rxy, I.gxy, I.Gxy and I.D) 
    % and (I.Gy, I.gy)
    block.I.gxy = block.X(p+1:T,:)'*Y(p+1:T,:);
    for n=1:N
        block.I.rxy(:,:,n) = dy(:,:,n)*X(p+1:T,:);
        for ki=1:k
            if block.p > 1
                Dtmp = dy(:,:,n)*squeeze(block.dX(:,ki,:))';
            else
                % With p=1, 'squeeze' already tranposes singleton dimension
                Dtmp = dy(:,:,n)*squeeze(block.dX(:,ki,:));
            end
            Dv = Dtmp(:)';
            block.I.D(ki,:,n) = Dv;
            if block.p > 1
                block.I.Gxy(:,ki,n) = squeeze(block.dX(:,ki,:))*Y(p+1:T,n);
            else
                block.I.Gxy(:,ki,n) = squeeze(block.dX(:,ki,:))'*Y(p+1:T,n);
            end
        end
        block.I.Gy(:,:,n) = dy(:,:,n)*dy(:,:,n)';
        block.I.gy(:,n)   = dy(:,:,n)*Y(p+1:T,n);
    end
end

%-Setting up spatial permutation matrices
%--------------------------------------------------------------------------
if block.verbose
    disp('Setting up spatial permutation matrices');
end
% Set up permutation matrix for regression coefficients
Nk = N*k;
block.Hw = sparse(Nk,Nk);
ii = [];
for kk=1:k
    ii = [ii, kk:k:Nk];
end
for nk=1:Nk
    block.Hw(ii(nk),nk) = 1;
end
% Set up permutation matrix for AR coefficients
if p > 0
    Np = N*p;
    block.Ha = sparse(Np,Np);
    ii = [];
    for pp=1:p
        ii = [ii, pp:p:Np];
    end
    for np=1:Np
        block.Ha(ii(np),np) = 1;
    end
end

try
    block.compute_det_D;
catch
    block.compute_det_D = 0;
end

if block.update_F
    if block.compute_det_D
        if block.verbose
            disp('Computing log determinant of spatial precision matrices');
        end
        % Get log determinant of Dw
        [vvv,ddd] = eig(full(block.Dw));
        dd = diag(ddd);
        dd = dd(dd>eps);
        block.log_det_Dw = sum(log(dd));
        if p > 0
            if isfield(block,'Da')
                if block.Da==block.Dw
                    block.log_det_Da = block.log_det_Dw;
                else
                    % Get log determinant of Da
                    [vvv,ddd] = eig(full(block.Da));
                    dd        = diag(ddd);
                    dd        = dd(dd>eps);
                    block.log_det_Da = sum(log(dd));
                end
            end
        end
    else
        block.log_det_Dw = 0;
        if p > 0
            block.log_det_Da = 0;
        end
    end
end

%-Computing data projections
%--------------------------------------------------------------------------
if block.verbose
    disp('Computing data projections');
end
% Set up design and data projections
try 
    block.XT;
    block.XTX;
catch
    block.XT  = X';
    block.XTX = block.XT*X;
end

for n=1:N
    block.XTY(:,n) = block.XT*Y(:,n);
    block.y2(n)    = Y(p+1:T,n)'*Y(p+1:T,n);
end

% Precompute quantities for the Negative Free Energy
block.C2 = N*T*log(2*pi);
