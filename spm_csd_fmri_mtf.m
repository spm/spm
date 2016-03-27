function [y,w,S,Gu,Gn] = spm_csd_fmri_mtf(P,M,U)
% Spectral response of a DCM (transfer function x noise spectrum)
% FORMAT [y,w,S,Gu,Gn] = spm_csd_fmri_mtf(P,M,U)
%
% P - model parameters
% M - model structure
% U - model inputs (expects U.csd as complex cross spectra)
%
% y - y(nw,nn,nn} - cross-spectral density for nn nodes
%                 - for nw frequencies in M.Hz
% w  - frequencies
% S  - directed transfer functions (complex)
% Gu - CSD of neuronal fluctuations
% Gn - CSD of observation noise
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
% $Id: spm_csd_fmri_mtf.m 6759 2016-03-27 19:45:17Z karl $


% compute log-spectral density
%==========================================================================

% frequencies of interest
%--------------------------------------------------------------------------
w    = M.Hz(:);
nw   = length(w);

% number of nodes and endogenous (neuronal) fluctuations
%--------------------------------------------------------------------------
nn   = size(P.A,1);
nu   = length(M.u);
form = '1/f';


% spectrum of neuronal fluctuations (Gu) and observation noise (Gn)
%==========================================================================

% experimental inputs
%--------------------------------------------------------------------------
Gu    = zeros(nw,nu,nu);
Gn    = zeros(nw,nn,nn);
if any(any(P.C))
    for i = 1:nu
        for j = 1:nu
            for k = 1:nw
                Gu(k,i,j) = P.C(i,:)*squeeze(U.csd(k,:,:))*P.C(j,:)';
            end
        end
    end
end

% neuronal fluctuations (Gu) (1/f or AR(1) form)
%--------------------------------------------------------------------------
for i = 1:nu
    if strcmp(form,'1/f')
        G     = w.^(-exp(P.a(2,1)));
    else
        G     = spm_mar2csd(exp(P.a(2,1))/2,w);
    end
    Gu(:,i,i) = Gu(:,i,i) + exp(P.a(1,1))*G/sum(G);
end

% region specific observation noise (1/f or AR(1) form)
%--------------------------------------------------------------------------
for i = 1:nn
    if strcmp(form,'1/f')
        G     = w.^(-exp(P.c(2,i))/2);
    else
        G     = spm_mar2csd(exp(P.c(2,i))/2,w);
    end
    Gn(:,i,i) = Gn(:,i,i) + exp(P.c(1,i))*G/sum(G);
end


% global components
%--------------------------------------------------------------------------
if strcmp(form,'1/f')
    G = w.^(-exp(P.b(2,1))/2);
else
    G = spm_mar2csd(exp(P.b(2,1))/2,w);
end
for i = 1:nn
    for j = i:nn
        Gn(:,i,j) = Gn(:,i,j) + exp(P.b(1,1))*G/sum(G);
        Gn(:,j,i) = Gn(:,i,j);
    end
end


% transfer functions (FFT of first-order Volterra kernel)
%==========================================================================
P.C   = speye(nn,nu);
S     = spm_dcm_mtf(P,M);

% predicted cross-spectral density
%--------------------------------------------------------------------------
G     = zeros(nw,nn,nn);
for i = 1:nw
    G(i,:,:) = reshape(S(i,:,:),nn,nn)*reshape(Gu(i,:,:),nn,nn)*reshape(S(i,:,:),nn,nn)';
end


% and channel noise
%--------------------------------------------------------------------------
if isfield(M,'g')
    y = G + Gn;
else
    y = G;
end

