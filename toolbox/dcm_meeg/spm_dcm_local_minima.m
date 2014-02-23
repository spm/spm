function spm_dcm_local_minima(DCM)
% evaluates the free energy landscape around the posterior
% FORMAT: spm_dcm_local_minima(DCM)
% DCM - (invert) model structure
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_local_minima.m 5892 2014-02-23 11:00:16Z karl $

% find dimension of greatest curvature
%==========================================================================
U  = spm_svd(DCM.Cp);
U  = U(:,end);
Cu = U'*DCM.Cp*U;
Su = spm_sqrtm(Cu);

% Free energy landscape
%==========================================================================
DCM.options.DATA = 0;
DCM.M.Nmax = 1;
DCM.name   = 'test';

N = 128;
s = linspace(-8,8,N);
for i = 1:N
    
    Ep           = spm_vec(DCM.Ep);
    Ep           = Ep + s(i)*U*Su*(U'*Ep);
    DCM.M.P      = spm_unvec(Ep,DCM.Ep);
    try
        [Qp,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.xU,DCM.xY);
        FF(i)        = F;
    catch
        FF(i)    = NaN;
    end
end

plot(s,FF)
title('Free energy','FontSize',16);
xlabel('parameter 1')
ylabel('parameter 2')
axis square





