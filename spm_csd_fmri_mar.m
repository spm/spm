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
% $Id: spm_csd_fmri_mar.m 5606 2013-08-13 08:47:22Z karl $


% compute log-spectral density
%==========================================================================

% lags of interest
%--------------------------------------------------------------------------
try, pst = M.pst(:);      end
try, pst = M.dt*(1:M.N)'; end

% number of nodes and endogenous (neuronal) fluctuations
%--------------------------------------------------------------------------
np   = M.p;                                   % number of MAR lags
nn   = M.l;                                   % number of nodes (regions)
nu   = length(M.u);                           % number of input (regions)


% cross-covaraince functions of neuronal fluctations (Vu) and noise (Vn)
%==========================================================================

% experimental inputs
%--------------------------------------------------------------------------
for i = 1:nu
    for j = 1:nu
        if size(P.C,2)
            for k = 1:nw
                Vu(i,j) = P.C(i,:)*squeeze(U.csd(k,:,:))*P.C(j,:)';
            end
        else
            Vu(i,j) = sparse(M.N,M.N);
        end
    end
end


% neuronal inputs
%--------------------------------------------------------------------------
for i = 1:nu
    Vu{i,i} = Vu{i,i} + exp(P.a(1,i))*spm_Q(exp(P.a(2,i))/2,M.N);
end
Vu    = spm_cat(Vu);


% observation noise
%--------------------------------------------------------------------------
for i = 1:nn
    
    % global component
    %----------------------------------------------------------------------
    for j = 1:nn
        V       = exp(P.b(1,1))*spm_Q(exp(P.b(2,1))/2,np + 1)/32;
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
P.C   = speye(nn,nu);
[S,k] = spm_dcm_mtf(P,M);

% matix form
%--------------------------------------------------------------------------
for i = 1:size(k,2)
    for j = 1:size(k,3)
        K{i,j} = k(:,i,j);
    end
end
K     = spm_cat(K);

% lagged matix form (with a decimation factor of R)
%--------------------------------------------------------------------------
R     = 4;
for p = 1:np
    t = (1 + p*R):size(k,1);
    for i = 1:size(k,2)
        for j = 1:size(k,3)
            L{i,j}    = zeros(size(k,1),1);
            L{i,j}(t) = k(t - p*R,i,j);
        end
    end
    KK{1,p} = spm_cat(L);
end
KK    = spm_cat(KK);

% predicted MAR coefficients
%--------------------------------------------------------------------------
A     = spm_inv(KK'*Vu*KK + Vn)*(KK'*Vu*K + Rn);




% predicted MAR coeficients
%--------------------------------------------------------------------------
G     = zeros(nw,nn,nn);
for i = 1:nn
    for j = 1:nn
        for k = 1:nu
            for l = 1:nu
                G(:,i,j) = G(:,i,j) + S(:,i,k).*Gu(:,k,l).*conj(S(:,j,l));
            end
        end
    end
end

% and  channel noise
%--------------------------------------------------------------------------
if isfield(M,'g')
    y = G + Gn;
else
    y = G;
end
