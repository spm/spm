function [slice] = spm_vb_a (Y,slice)
% Update AR coefficients in VB GLM-AR model 
% FORMAT [slice] = spm_vb_a (Y,slice)
%
% Y             [T x N] time series 
% slice         data structure 
%
% @(#)spm_vb_a.m	2.1 Will Penny and Nelson Trujillo-Barreto 04/08/04

if slice.verbose
    disp('Updating a');
end

p=slice.p;
k=slice.k;
N=slice.N;
T=slice.T;
Jk = kron(diag(slice.mean_beta),slice.D);
J=slice.Ha*Jk*slice.Ha';

Y_err_w = Y(p+1:T,:) - slice.X(p+1:T,:)*reshape(slice.w_mean,k,N);

for n=1:N,
    
    % Set up indexes
    wblock_n      = [(n-1)*k+1:n*k]; % Indexes for regression coeffs
    block_n       = [(n-1)*p+1:n*p]; % Indexes for AR coeffs
    block_ni          = [1:N*p];
    block_ni(block_n) = [];
    Jnn               = J(block_n,block_n);
    Jni               = J(block_n,block_ni);
    
    w_mean=slice.w_mean(wblock_n);
    
    % C_tilde and D_tilde:
    % Previously equation 50 in paper VB1
    % but now using more efficient cross-covariance formulae in paper VB3
    C1=slice.I.Gy(:,:,n);
    C2=reshape(slice.I.S'*slice.w2{n}(:),p,p);
    C3=-slice.I.W_tilde(:,:,n);
    C4=C3';
    C_til=C1+C2+C3+C4;
   
    D1=slice.I.gy(:,n)';
    D2=(-slice.I.rxy(:,:,n)*w_mean)';
    D3=(-slice.I.Gxy(:,:,n)*w_mean)';
    D4=slice.w2{n}(:)'*slice.I.R1;
    D_til=D1+D2+D3+D4;
    
    % Spatially regularised AR coefficients - paper VB4
    slice.a_cov{n}  = inv(slice.mean_lambda(n)*C_til + Jnn);
    a_mean = slice.a_cov{n}*(slice.mean_lambda(n)*D_til'-Jni*slice.a_mean(block_ni,1));
    slice.a_mean(block_n,1)=a_mean;
    slice.a2{n}=a_mean*a_mean'+slice.a_cov{n};
    
    % Intermediate quantities for updating other parameters
    % - see section on cross-covariance formulae in paper VB3
    slice.I.A2_tilde(:,:,n)=reshape(slice.I.S*slice.a2{n}(:),k,k);
    slice.I.A3a_tilde(:,:,n)=-reshape(slice.I.R1*a_mean,k,k);
end
 
slice.ap_mean = reshape(slice.a_mean,p,N);

