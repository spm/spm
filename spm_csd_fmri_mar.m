function [y,S,k] = spm_csd_fmri_mar(P,M,U)
% Prediction of MAR coefficients for DCM
% FORMAT [y,S,K] = spm_csd_fmri_mar(P,M,U)
%
% P - model parameters
% M - model structure
% U - model inputs (expects U.csd as complex cross spectra)
%
% y - y(nw,nn,nn} - cross-spectral density for nn nodes
%                 - for nw frequencies in M.Hz
% K - Volterra kernels
% S - directed transfer functions (complex)
%
% This routine computes the spectral response of a network of regions
% driven by  endogenous fluctuations and exogenous (experimental) inputs.
% It returns the complex cross spectra of regional responses as a
% three-dimensional array. The endogenous innovations or fluctuations are
% parameterised in terms of a (scale free) power law, in frequency space.
%
% When the observer function M.g is specified, the CSD response is
% supplemented with observation noise in sensor space; otherwise the CSD
% is noiseless.
%
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_csd_fmri_mar.m 5618 2013-08-17 10:36:56Z karl $


% number of nodes and endogenous (neuronal) fluctuations
%--------------------------------------------------------------------------
np   = M.p;                                   % number of MAR lags
nn   = M.l;                                   % number of nodes (regions)

% cross-covaraince functions of neuronal fluctations (Vu) and noise (Vn)
%==========================================================================

% experimental inputs
%--------------------------------------------------------------------------
for i = 1:nn
    for j = 1:nn
        if any(any(P.C))
            for k = 1:M.N
                V(k) = P.C(i,:)*squeeze(U.ccf(k,:,:))*P.C(j,:)';
            end
            Vu{i,j} = toeplitz(V);
        else
            Vu{i,j} = sparse(M.N,M.N);
        end
    end
end


% neuronal inputs
%--------------------------------------------------------------------------
for i = 1:nn
    Vu{i,i} = Vu{i,i} + exp(P.a(1,i))*spm_Q(exp(P.a(2,i))/2,M.N);
end
Vu    = spm_cat(Vu);


% observation noise
%--------------------------------------------------------------------------
for i = 1:nn
    
    % global component
    %----------------------------------------------------------------------
    for j = 1:nn
        V       = exp(P.b(1,1))*spm_Q(exp(P.b(2,1))/2,np + 1)/64;
        Vn{i,j} = V((1:np),(1:np));
        Rn{i,j} = V((1:np) + 1,1);
    end
    
    % region specific
    %----------------------------------------------------------------------
    V       = exp(P.c(1,i))*spm_Q(exp(P.c(2,i))/2,np + 1)/8;
    Vn{i,i} = Vn{i,j} + V((1:np),(1:np));
    Rn{i,i} = Rn{i,j} + V((1:np) + 1,1);
    
end
Vn    = spm_cat(Vn);
Rn    = spm_cat(Rn);


% first-order Volterra kernel
%==========================================================================
P.C   = speye(nn,nn);
[S,k] = spm_dcm_mtf(P,M);

% matix form
%--------------------------------------------------------------------------
for i = 1:nn
    for j = 1:nn
        K{i,j} = k(:,i,j);
    end
end
K     = spm_cat(K);

% lagged matix form
%--------------------------------------------------------------------------
for i = 1:nn
    for j = 1:nn
        KK{i,j}    = zeros(nn,np);
        for p = 1:np
            t = (1 + p):size(k,1);
            KK{i,j}(t,p) = k(t - p,i,j);
        end
    end
end
KK    = spm_cat(KK);

% predicted MAR coefficients
%--------------------------------------------------------------------------
y     = spm_inv(KK'*Vu*KK + Vn)*(KK'*Vu*K + Rn);
