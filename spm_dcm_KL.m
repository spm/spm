function [D,C,K] = spm_dcm_KL(M)
% Computes the distance between two models based on prior responses
% FORMAT [D,C,K] = spm_dcm_KL(Mi,Mj)
%
% M{1:n}   - structure array of models
%
% D(n x n) - distance matrix (KL divergence)
% C{1:n}   - response covariances
% K{1:n}   - response means
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_KL.m 2804 2009-03-02 12:03:00Z karl $
 
 
% Volterra kernels
%==========================================================================
 
% time bins (if not specified)
%--------------------------------------------------------------------------
try 
    dt = M{1}.dt;
    N  = M{1}.N;
catch
    
    % Bilinear representation
    %----------------------------------------------------------------------
    M0 = spm_bireduce(M{1},P{1});
    s  = real(eig(full(M0)));
    s  = max(s(find(s < 0)));
    N  = 32;
    dt = -4/(s*N);
    
end
 
% get covariances of prior responses
%==========================================================================
m     = length(M);
for i = 1:m
    
    % Get parameters (adding a little to prevent expansion around zero)
    %----------------------------------------------------------------------
    P      = M{i}.pE;
    pC     = M{i}.pC;
    P      = spm_vec(P) + sqrt(diag(pC))/8;
    P      = spm_unvec(P,M{i}.pE);
    
    % add a little to inputs (to prevent expansion around zero)
    %----------------------------------------------------------------------
    M{i}.u = ones(M{i}.m,1)/8;
    
    % get eigen-space of parameters for computational efficiency
    %----------------------------------------------------------------------
    V      = spm_svd(pC);
    pC     = V'*pC*V;
    
    % get derivative of kernels w.r.t. parameters
    %----------------------------------------------------------------------
    [dkdp k] = spm_diff('spm_kernel',M{i},P,N,dt,2,{[],V});
    
    % prior mean and covariance of kernels
    %----------------------------------------------------------------------
    k     = spm_vec(k);
    dk    = sparse(length(k),length(dkdp));
    for j = 1:length(dkdp)
        dk(:,j) = spm_vec(dkdp{j});
    end
    K{i} = k;
    C{i} = dk*pC*dk';
    n(i) = length(pC);
    
end
 
 
% evaluate KL divergence
%--------------------------------------------------------------------------
for i = 1:m
    for j = 1:m
        Pi     = spm_pinv(C{i});
        k      = K{i} - K{j};
        d      = spm_logdet(C{i}) - spm_logdet(C{j}) + ...
                 trace(Pi*C{j}) + k'*Pi*k - n(i);
        D(i,j) = d/2;
    end
end

% ensure D is symmetric
%--------------------------------------------------------------------------
D   = (D + D')/2;

