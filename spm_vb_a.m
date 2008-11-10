function [block] = spm_vb_a (Y,block)
% Update AR coefficients in VB GLM-AR model 
% FORMAT [block] = spm_vb_a (Y,block)
%
% Y             [T x N] time series 
% block         data structure 
%
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny and Nelson Trujillo-Barreto
% $Id: spm_vb_a.m 2451 2008-11-10 16:20:32Z lee $

if block.verbose
    disp('Updating a');
end

p=block.p;
k=block.k;
N=block.N;
T=block.T;

if strcmp(block.priors.A,'Discrete')
    as=zeros(block.p,block.priors.S);
else
    Jk = kron(diag(block.mean_beta),block.Da);
    J=block.Ha*Jk*block.Ha';
end

Y_err_w = Y(p+1:T,:) - block.X(p+1:T,:)*reshape(block.w_mean,k,N);

for n=1:N,
    
    % Set up indexes
    subblock_n      = [(n-1)*k+1:n*k]; % Indexes for regression coeffs
    w_mean=block.w_mean(subblock_n);
    
    % C_tilde and D_tilde:
    % Previously equation 50 in paper VB1
    % but now using more efficient cross-covariance formulae in paper VB3
    C1=block.I.Gy(:,:,n);
    C2=reshape(block.I.S'*block.w2{n}(:),p,p);
    C3=-block.I.W_tilde(:,:,n);
    C4=C3';
    C_til=C1+C2+C3+C4;
   
    D1=block.I.gy(:,n)';
    D2=(-block.I.rxy(:,:,n)*w_mean)';
    D3=(-block.I.Gxy(:,:,n)*w_mean)';
    D4=block.w2{n}(:)'*block.I.R1;
    D_til=D1+D2+D3+D4;
    
    block_n       = [(n-1)*p+1:n*p]; % Indexes for AR coeffs
    % Spatially regularised AR coefficients - paper VB4
    switch block.priors.A
        case 'Discrete',
            beta_n=diag(block.priors.gamma(n,:)*block.mean_beta');
            a_n=block.priors.gamma(n,:)*block.as';
            block.a_cov{n}  = inv(block.mean_lambda(n)*C_til + beta_n);
            a_mean = block.a_cov{n}*(block.mean_lambda(n)*D_til'+beta_n*a_n');
            as = as + a_mean*block.priors.gamma(n,:);
        otherwise,
            block_ni          = [1:N*p];
            block_ni(block_n) = [];
            Jnn               = J(block_n,block_n);
            Jni               = J(block_n,block_ni);
            block.a_cov{n}  = inv(block.mean_lambda(n)*C_til + Jnn);
            a_mean = block.a_cov{n}*(block.mean_lambda(n)*D_til'-Jni*block.a_mean(block_ni,1));
    end
    block.a_mean(block_n,1)=a_mean;
    block.a2{n}=a_mean*a_mean'+block.a_cov{n};
    
    % Intermediate quantities for updating other parameters
    % - see section on cross-covariance formulae in paper VB3
    block.I.A2_tilde(:,:,n)=reshape(block.I.S*block.a2{n}(:),k,k);
    block.I.A3a_tilde(:,:,n)=-reshape(block.I.R1*a_mean,k,k);
end
 
block.ap_mean = reshape(block.a_mean,p,N);

if strcmp(block.priors.A,'Discrete')
    % Update AR means for each category
    for s=1:block.priors.S,
        block.as(:,s)=as(:,s)/block.priors.N(s);
    end
end

