function [Gu,Gs,Gn,f] = spm_csd_mtf_gu(P,M)
% Spectral desnities of innovations and noise for DCM for CSD
% FORMAT [Gu,Gs,Gn,f] = spm_csd_mtf_gu(P,M)
% FORMAT [Gu,Gs,Gn,f] = spm_csd_mtf_gu(P,f)
%
% P   - parameters
% M   - neural mass model structure (with M.Hz)
% f   - frequencies of interest (Hz)
%
% Gu  - neuronal innovations
% Gn  - channel noise (non-specific)
% Gs  - channel noise (specific)
%
% f   - frequency
%
% fluctuations and noise parameters: for n regions and c channels
%--------------------------------------------------------------------------
%    pE.a(2,n) - neuronal fluctuations        - amplitude and exponent
%    pE.b(2,c) - channel noise (non-specific) - amplitude and exponent
%    pE.c(2,c) - channel noise (specific)     - amplitude and exponent
%    pE.d(8,n) - neuronal fluctuations        - basis set coefficients
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_csd_mtf_gu.m 6856 2016-08-10 17:55:05Z karl $

 
% frequencies of interest
%--------------------------------------------------------------------------
try, f = M.Hz(:); catch, f = M(:); end

% Number of sources and fequencies
%--------------------------------------------------------------------------
if isfield(P,'d')
    ns = max(size(P.a,2),size(P.d,2));
else
    ns = size(P.a,2);
end
nf     = size(f,1);


% spectrum of neuronal innovations (Gu)
%==========================================================================
for i = 1:size(P.a,2)
    Gu(:,i) = exp(P.a(1,i))*f.^(-exp(P.a(2,i)));
end

if size(Gu,2) == 1, Gu = Gu*ones(1,ns); end


% add structured innovations - a discrete cosine set of order length(P.d)
%--------------------------------------------------------------------------
if isfield(P,'d')
    nd = size(P.d,1);
    X  = spm_dctmtx(nf,nd + 1);
    Mu = exp(X(:,2:end)*P.d);
else
    Mu = ones(nf,1);
end

if size(Mu,2) == 1, Mu = Mu*ones(1,ns); end

% return unless channel noise is required
%--------------------------------------------------------------------------
Gu     = Gu.*Mu;
if nargout < 2, return, end


% spectrum of channel noise (non-specific)
%==========================================================================
Gn  = exp(P.b(1) - 2)*f.^(-exp(P.b(2))); 

% and spectrum of channel noise (specific: with the same exponent)
%--------------------------------------------------------------------------
for i = 1:size(P.c,2)
    Gs(:,i) = exp(P.c(1,i) - 2)*f.^(-exp(P.c(2,1)));
end



