function [slice] = spm_vb_init_slice (Y,slice)
% Initialise Variational Bayes for GLM-AR models
% FORMAT [slice] = spm_vb_init_slice (Y,slice)
%
% Y             [T x N] time series with T time points, N voxels
% slice         GLM-AR data structure
%
% %W% Will Penny and Nelson Trujillo-Barreto %E%

k=slice.k;
p=slice.p;
N=slice.N;
T=slice.T;
X=slice.X;

%% Default optimisation parameters
if ~isfield(slice,'tol')
    slice.tol=0.0001;
end
if ~isfield(slice,'maxits')
    slice.maxits=4;
end
if ~isfield(slice,'verbose')
    slice.verbose=1;
end

if slice.verbose
    disp('Initialising slice');
    disp(' ');
end

%% Default priors
if ~isfield(slice,'b_alpha_prior')
    slice.b_alpha_prior=10;
end   
if ~isfield(slice,'c_alpha_prior')
    slice.c_alpha_prior=0.1;
end
if ~isfield(slice,'b_beta_prior')
    slice.b_beta_prior=10;
end   
if ~isfield(slice,'c_beta_prior')
    slice.c_beta_prior=0.1;
end

if ~isfield(slice,'b_lambda_prior')
    slice.b_lambda_prior=10;
end
if ~isfield(slice,'c_lambda_prior')
    slice.c_lambda_prior=0.1;
end


% Initialise approximate alpha posterior
if isfield(slice,'mean_alpha')
    k1=size(slice.mean_alpha);
    if ~(k1==k)
        disp('Error in spm_vbglmar_slice: mean_alpha incorrect dimension');
        return
    end
    slice.b_alpha=slice.mean_alpha/slice.c_alpha_prior;
    slice.c_alpha=slice.c_alpha_prior*ones(k,1);
else
    slice.b_alpha=slice.b_alpha_prior*ones(k,1);
    slice.c_alpha=slice.c_alpha_prior*ones(k,1);
    slice.mean_alpha=slice.b_alpha.*slice.c_alpha;
end
slice.c_alpha  = (N./2 + slice.c_alpha_prior)*ones(k,1);

% Initialise approximate beta posterior
if isfield(slice,'mean_beta')
    p1=size(slice.mean_beta);
    if ~(p1==p)
        disp('Error in spm_vbglmar_slice: mean_beta incorrect dimension');
        return
    end
    slice.b_beta=slice.mean_beta/slice.c_beta_prior;
    slice.c_beta=slice.c_beta_prior*ones(p,1);
else
    slice.b_beta=slice.b_beta_prior*ones(p,1);
    slice.c_beta=slice.c_beta_prior*ones(p,1);
    slice.mean_beta=slice.b_beta.*slice.c_beta;
end
slice.c_beta  = (p/2 + slice.c_beta_prior)*ones(p,1);

% Initialise approximate lambda posterior
if isfield(slice,'mean_lambda')
    n1=size(slice.mean_lambda,1);
    if ~(n1==N)
        disp('Error in spm_vbglmar_slice: mean_lambda incorrect dimension');
        return
    end
    slice.b_lambda=slice.mean_lambda/slice.c_lambda_prior;
    slice.c_lambda=slice.c_lambda_prior*ones(N,1);
else
    slice.b_lambda=slice.b_lambda_prior*ones(N,1);
    slice.c_lambda=slice.c_lambda_prior*ones(N,1);
    slice.mean_lambda=slice.b_lambda.*slice.c_lambda;
end
slice.c_lambda = ((T-slice.p)./2 + slice.c_lambda_prior)*ones(N,1);

disp('Initialising regression coefficient posterior');
% Initialise approximate w posterior
try
    Xp=slice.Xp;
    X2=slice.X2;
catch
    [ux,dx,vx]=svd(X);
    ddx=diag(dx);
    svd_tol=max(ddx)*eps*k;
    rank_X=sum(ddx > svd_tol);
    ddxm=diag(ones(rank_X,1)./ddx(1:rank_X));
    ddxm2=diag(ones(rank_X,1)./(ddx(1:rank_X).^2));
    Xp=vx(:,1:rank_X)*ddxm*ux(:,1:rank_X)';
    X2=vx(:,1:rank_X)*ddxm2*vx(:,1:rank_X)';
end  

w_ols=Xp*Y;
w_mean = w_ols;
Y_pred = X*w_mean;
v=mean((Y-Y_pred).^2);
for n=1:N,
    w_cov_temp = v(n)*X2;
    slice.w_cov{n} = w_cov_temp;
end;
slice.w_mean=w_mean(:);
slice.w_ols=slice.w_mean;
slice.wk_mean      = reshape(slice.w_mean,k,N);
slice.wk_ols      = reshape(slice.w_ols,k,N);

% Initialise AR coefficient posterior
disp('Initialising AR coefficients');
% Embed data
for pp=1:p,
    dy(pp,1:T-p,:)=Y(p-pp+1:T-pp,:);
end
if p>0
    e = Y(p+1:T,:) - Y_pred(p+1:T,:);
    for n=1:N,
        for pp=1:p,
            dyhat(pp,:)=w_mean(:,n)'*squeeze(slice.dX(pp,:,:));
        end
        E_tilde=dy(:,:,n)-dyhat;
        iterm       = inv(E_tilde * E_tilde');
        slice.ap_mean(:,n) = (iterm * E_tilde*e(:,n));  
        e_pred      = E_tilde' * slice.ap_mean(:,n);
        v2          = mean((e(:,n) - e_pred).^2);
        slice.a_cov{n} = v2 * iterm;
        slice.a2{n}=slice.ap_mean(:,n)*slice.ap_mean(:,n)'+slice.a_cov{n};
    end
    slice.ap_ols=slice.ap_mean;
    slice.a_mean=slice.ap_mean(:);
end

if p>0
    disp('Setting up cross-covariance matrices');
    % Get input-output lagged covariances (I.rxy, I.gxy, I.Gxy and I.D) 
    % and (I.Gy, I.gy)
    slice.I.gxy=slice.X(p+1:T,:)'*Y(p+1:T,:);
    for n=1:N,
        slice.I.rxy(:,:,n)=dy(:,:,n)*X(p+1:T,:);
        for ki=1:k,
            if slice.p>1
                Dtmp=dy(:,:,n)*squeeze(slice.dX(:,ki,:))';
            else
                % With p=1, 'squeeze' already tranposes singleton dimension
                Dtmp=dy(:,:,n)*squeeze(slice.dX(:,ki,:));
            end
            Dv=Dtmp(:)';
            slice.I.D(ki,:,n)=Dv;
            if slice.p>1
                slice.I.Gxy(:,ki,n)=squeeze(slice.dX(:,ki,:))*Y(p+1:T,n);
            else
                slice.I.Gxy(:,ki,n)=squeeze(slice.dX(:,ki,:))'*Y(p+1:T,n);
            end
        end
        slice.I.Gy(:,:,n)=dy(:,:,n)*dy(:,:,n)';
        slice.I.gy(:,n)=dy(:,:,n)*Y(p+1:T,n);
    end
end

disp('Setting up spatial permutation matrices');
% Set up permutation matrix for regression coefficients
Nk = N*k;
slice.Hw=sparse(Nk,Nk);
ii=[];
for kk=1:k,
    ii=[ii, kk:k:Nk];
end
for nk=1:Nk,
    slice.Hw(ii(nk),nk)=1;
end
% Set up permutation matrix for AR coefficients
if p > 0
    Np = N*p;
    slice.Ha=sparse(Np,Np);
    ii=[];
    for pp=1:p,
        ii=[ii, pp:p:Np];
    end
    for np=1:Np,
        slice.Ha(ii(np),np)=1;
    end
end

if slice.update_F
    disp('Computing log determinant of spatial precision matrix for evidence');
    % Get log determinant of D
    [vvv,ddd]=eig(full(slice.D));
    slice.log_det_D=sum(log(diag(ddd)));
end

disp('Computing data projections');
% Set up design and data projections
try 
    slice.XT;
    slice.XTX;
catch
    slice.XT=X';
    slice.XTX=slice.XT*X;
end

for n=1:N,
    slice.XTY(:,n)=slice.XT*Y(:,n);
    slice.y2(n)=Y(p+1:T,n)'*Y(p+1:T,n);
end

% Precompute quantities for the Negative Free Energy
slice.C2  = N*T*log(2*pi);


