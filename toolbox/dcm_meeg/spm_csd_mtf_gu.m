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
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_csd_mtf_gu.m 4814 2012-07-30 19:56:05Z karl $

 
% frequencies of interest
%--------------------------------------------------------------------------
try, f = M.Hz(:); catch, f = M(:); end


% spectrum of neuronal innovations (Gu)
%==========================================================================
for i = 1:size(P.a,2)
    Gu(:,i) = exp(P.a(1,i))*f.^(-exp(P.a(2,i)));
end

% add structured innovations - a discrete cosine set of order length(P.d)
%--------------------------------------------------------------------------
if isfield(P,'d')
    X  = spm_dctmtx(size(Gu,1),size(P.d,1) + 1);
    Gu = Gu.*exp(X(:,2:end)*P.d);
end

if size(Gu,2) == 1, Gu = Gu*ones(1,size(P.a,2)); end

% return unless channel noise is required
%--------------------------------------------------------------------------
if nargout < 2, return, end


% spectrum of channel noise (non-specific)
%==========================================================================
Gn  = exp(P.b(1))*f.^(-exp(P.b(2)))*4; 

% and spectrum of channel noise (specific)
%--------------------------------------------------------------------------
for i = 1:size(P.c,2)
    Gs(:,i) = exp(P.c(1,i))*f.^(-exp(P.c(2,i)))*4;
end



