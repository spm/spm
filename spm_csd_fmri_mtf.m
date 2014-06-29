function [y,w,S] = spm_csd_fmri_mtf(P,M,U)
% Spectral response of a DCM (transfer function x noise spectrum)
% FORMAT [y,w,s] = spm_csd_fmri_mtf(P,M,U)
%
% P - model parameters
% M - model structure
% U - model inputs (expects U.csd as complex cross spectra)
%
% y - y(nw,nn,nn} - cross-spectral density for nn nodes
%                 - for nw frequencies in M.Hz
% w - frequencies
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
% $Id: spm_csd_fmri_mtf.m 6075 2014-06-29 21:11:40Z karl $


% compute log-spectral density
%==========================================================================

% frequencies of interest
%--------------------------------------------------------------------------
w    = M.Hz(:);
nw   = length(w);

% number of nodes and endogenous (neuronal) fluctuations
%--------------------------------------------------------------------------
nn   = M.l;
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

% amplitude scaling constant
%--------------------------------------------------------------------------
C     = 1/256;

% neuronal fluctuations (Gu) (1/f or AR(1) form)
%--------------------------------------------------------------------------
for i = 1:nu
    if strcmp(form,'1/f')
        G     = exp(P.a(1,i))*w.^(-exp(P.a(2,i)))*4;
    else
        G     = exp(P.a(1,i))*spm_mar2csd(exp(P.a(2,i))/2,w,M.ns);
    end
    Gu(:,i,i) = Gu(:,i,i) + C*G;
end

% observation noise (1/f or AR(1) form)
%--------------------------------------------------------------------------
for i = 1:nn
    
    % global component
    %----------------------------------------------------------------------
    for j = 1:nn
        
        if strcmp(form,'1/f')
            G     = exp(P.b(1,1))*w.^(-exp(P.b(2,1))/2)/8;
        else
            G     = exp(P.b(1,1))*spm_mar2csd(exp(P.b(2,1))/2,w,M.ns)/64;
        end
        Gn(:,i,j) = Gn(:,i,j) + C*G;
    end
    
    % region specific
    %----------------------------------------------------------------------
    if strcmp(form,'1/f')
        G     = exp(P.c(1,i))*w.^(-exp(P.c(2,i))/2);
    else
        G     = exp(P.c(1,i))*spm_mar2csd(exp(P.c(2,i))/2,w,M.ns)/8;
    end
    Gn(:,i,i) = Gn(:,i,i) + C*G;
    
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

