function [pE,pC] = spm_hdm_priors(m,h)
% returns priors for a hemodynamic dynamic causal model
% FORMAT [pE,pC] = spm_hdm_priors(m,[h])
% m   - number of inputs
% h   - number of hemodynamic modes (default = 3)
%
% pE  - prior expectations
% pC  - prior covariances
%
% (5) biophysical parameters
%    P(1) - signal decay                  d(ds/dt)/ds)
%    P(2) - autoregulation                d(ds/dt)/df)
%    P(3) - transit time                  (t0)
%    P(4) - exponent for Fout(v)          (alpha)
%    P(5) - resting oxygen extraction     (E0)
%    P(6) - ratio of intra- to extra-     (epsilon)
%           vascular components of the
%           gradient echo signal   
%
% plus (m) efficacy priors
%    P(7) - ....
%
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_hdm_priors.m 4579 2011-12-02 20:21:07Z karl $



% default: 3 hemodynamic [eigen]modes
%---------------------------------------------------------------------------
if nargin < 2
    h = 3;
end

% biophysical parameters with prior expectation and
%---------------------------------------------------------------------------
pE    = [   0.65      0.41      0.98      0.32      0.34  ];

% covariance restricted to h modes (v) scaled by eigenvales (e) {see below)
%---------------------------------------------------------------------------
v     = [
        -0.0535    0.0095   -0.1117   -0.0040   -0.0026
        -0.0604   -0.0319    0.0430   -0.0077    0.0026
         0.1116   -0.0347   -0.2539   -0.0169   -0.0115
         0.1985    0.1698    0.4984   -0.4493    0.4434
         0.0029    0.2081    1.9582   -0.5209   -1.1634]';

e     = [2.1225    1.2006    0.3519    0.0039    0.0012];

% set variance of minor modes to zero
%---------------------------------------------------------------------------
i     = (h + 1):5;
e(i)  = 0;
pC    = v*diag(e)*v'/32;

% append scaling parameter for epsilon: prior variance = 1/32,
%---------------------------------------------------------------------------
pE    = [pE -1];
pC    = blkdiag(pC,1/128);

% append m efficacy priors
%---------------------------------------------------------------------------
pE    = [pE(:); zeros(m,1)];
pC    = blkdiag(pC,eye(m,m)*32);

return


% NOTES: sample covariances from Friston et al (2000)
%---------------------------------------------------------------------------
qC    = [0.0150    0.0052    0.0283    0.0002   -0.0027
         0.0052    0.0020    0.0104    0.0004   -0.0013
         0.0283    0.0104    0.0568    0.0010   -0.0069
         0.0002    0.0004    0.0010    0.0013   -0.0010
        -0.0027   -0.0013   -0.0069   -0.0010    0.0024];


% NOTES: Reduce rank of prior covariances for computational expediancy
%---------------------------------------------------------------------------

% assume independent priors in parameter space
%---------------------------------------------------------------------------
qC    = diag(diag(qC));


% model specification (single node DCM)
%---------------------------------------------------------------------------
M.f   = 'spm_fx_hdm';
M.g   = 'spm_gx_hdm';
M.x   = [0 0 0 0]';
M.pE  = [pE(1:6) 1];
M.m   = 1;
M.n   = 4;
M.l   = 1;
M.N   = 32;
M.dt  = 1/2;

% compute partial derivatives w.r.t. hemodynamic parameters [J] dy(t)/dp
%---------------------------------------------------------------------------
P     = M.pE;
p     = length(P);
dp    = 1e-6;
[k J] = spm_nlsi(M);
for i = 1:5
    M.pE    = P;
    M.pE(i) = M.pE(i) + dp;
    [k q]   = spm_nlsi(M);
    Jq(:,i) = (q - J)/dp;
end

% implied covariance of impulse response
%---------------------------------------------------------------------------
Cq    = Jq*qC*Jq';

% reduce to h hemodynamic modes in measurement space
%---------------------------------------------------------------------------
[v e] = spm_svd(Cq);
e     = diag(e);
v     = pinv(Jq)*v;
qC    = v*diag(e)*v';


% NOTES: graphics - eigenvalues of qC
%---------------------------------------------------------------------------
subplot(2,2,1)
bar(e)
xlabel('eigen mode')
title('eigenvalue')
set(gca,'XLim',[0 6])
axis square
grid on

% graphics - response differentials
%---------------------------------------------------------------------------
subplot(2,2,2)
plot([1:M.N]*M.dt,Jq*v(:,1),[1:M.N]*M.dt,Jq*v(:,2),'-.')
xlabel('PST {secs}')
title('hemodynamic modes')
axis square
grid on
