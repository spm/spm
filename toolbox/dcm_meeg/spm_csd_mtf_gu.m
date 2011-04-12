function [Gu,Gs,Gn,f] = spm_csd_mtf_gu(P,M)
% Spectral desnities of innovations and noise for DCM for CSD
% FORMAT [Gu,Gs,Gn,f] = spm_csd_mtf_gu(P,M)
%
% P   - parameters
% M   - neural mass model structure
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
% $Id: spm_csd_mtf_gu.m 4305 2011-04-12 18:15:32Z karl $

 
% frequencies of interest
%--------------------------------------------------------------------------
try
    dt = 1/(2*round(M.Hz(end)));
    N  = 1/dt;
catch
    N  = 128;
end
f     = [1:N/2]';


% spectrum of innovations (Gu)
%--------------------------------------------------------------------------
for i = 1:size(P.a,2)
    Gu(:,i) = exp(P.a(1,i))*f.^(-exp(P.a(2,i)));
end

% spectrum of channel noise (non-specific)
%--------------------------------------------------------------------------
Gn  = exp(P.b(1))*f.^(-exp(P.b(2)))*4; 


% spectrum of channel noise (specific)
%--------------------------------------------------------------------------
for i = 1:size(P.c,2)
    Gs(:,i) = exp(P.c(1,i))*f.^(-exp(P.c(2,i)))*4;
end
 
% add structured innovations - a discrete cosine set of order length(P.d)
%--------------------------------------------------------------------------
try
    X  = spm_dctmtx(N/2,size(P.d,1) + 1);
    Gu = Gu.*exp(X(:,2:end)*P.d);
end
if size(Gu,2) == 1, Gu = Gu*ones(1,M.m); end
if size(Gs,2) == 1, Gs = Gs*ones(1,M.l); end


