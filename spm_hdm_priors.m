function [pE,pC] = spm_hdm_priors(m)
% returns priors for a hemodynamic dynmaic causal model
% FORMAT [pE,pC] = spm_hdm_priors(m)
% m.. - number of inputs
%
% pE  - prior expectations
% pC  - prior covariances
%
% (5) biophysical parameters
%    P(1) - signal decay      - d(ds/dt)/ds)
%    P(2) - autoregulation    - d(ds/dt)/df)
%    P(3) - transit time                (t0)
%    P(4) - exponent for Fout(v)     (alpha)
%    P(5) - resting oxygen extraction   (E0)
%
% plus (m) efficacy priors
%    P(6) - ....
%
%___________________________________________________________________________
% %W% Karl Friston %E%

% biophysical parameters with prior expectation and
%---------------------------------------------------------------------------
pE    = [0.65      0.41      0.98      0.32      0.34  ];

% covariance restricted to 2 modes (v) scaled by eigenvales (e) {see below)
%---------------------------------------------------------------------------
v     = [-0.3833    0.0979   -0.9178   -0.0261   -0.0228;
          0.8604    0.3809   -0.3210    0.1041   -0.0268]';

e     =  [0.0332    0;
          0    0.0076];

pC    = v*e*v';

% append m efficacy priors
%---------------------------------------------------------------------------
pE    = [pE(:); zeros(m,1)];
pC    = blkdiag(pC,eye(m));

return


% NOTES: sample covariances from Friston et al (2000)
%---------------------------------------------------------------------------
qC    = [   0.0150    0.0052    0.0283    0.0002   -0.0027
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
M.fx  = 'spm_fx_HRF';
M.lx  = 'spm_lambda_HRF';
M.x   = [0 1 1 1]';
M.pE  = [pE 1];
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
h     = 1:2;					
[v s] = spm_svd(Cq);
v     = v(:,h);
Cq    = v*v'*Cq*v*v';
qC    = pinv(Jq)*Cq*pinv(Jq');


% NOTES: graphics - eigenvalues of qC
%---------------------------------------------------------------------------
subplot(2,2,1)
bar(diag(s))
xlabel('eigen mode')
title('eigenvalue')
set(gca,'XLim',[0 6])
axis square
grid on

% graphics - response differentials
%---------------------------------------------------------------------------
subplot(2,2,2)
plot([1:M.N]*M.dt,v(:,1),[1:M.N]*M.dt,v(:,2),'-.')
xlabel('PST {secs}')
title('hemodynamic modes')
axis square
grid on
