function [pE,pC,qE,qC] = spm_dcm_priors(A,B,C)
% returns the priors for a hemodynamic dynmaic causal model
% FORMAT [pE,pC,qE,qC] = spm_dcm_priors(A,B,C)
% A,B,C. - constraints on connections (1 - present, 0 - absent)
%
% pE     - prior expectations (connections and hemodynamic)
% pC     - prior covariances  (connections and hemodynamic)
% qE     - prior expectations (hemodynamic)
% qC     - prior covariances  (hemodynamic)
%___________________________________________________________________________
% %W% Karl Friston %E%

% number of regions
%---------------------------------------------------------------------------
n     = length(A);

% CONNECTIITY PRIORS
%===========================================================================
% covariances            pC = 1/((n - 1)*q)  - if A == 1
%                        pC = 0                if A == 0
%
% where, for aij = aji = 1/(n - 1)  => max(eig(J(0))) = 0
% and q is the Chi-squared threshold with n*(n - 1) df

% log(2)/b = half-life {b = self inhibition}
%---------------------------------------------------------------------------
p     = 1e-3;
b     = 1;
s     = spm_invNcdf(1 - p);
q     = spm_invXcdf(1 - p,n*(n - 1));
q     = n/((n - 1)*q);

% intrinsic connections A {additional priors from eigenvalues}
%---------------------------------------------------------------------------
A     = A - diag(diag(A));
pC    = diag([(b/s)^2; A(:)*q; B(:)*q; C(:)]);

% expectations
%---------------------------------------------------------------------------
A     = -speye(n,n);
B     = B*0;
C     = C*0;
pE    = [b; A(:); B(:); C(:)];


% HEMODYNAMIC PRIORS
%===========================================================================
% P(1) - signal decay     - d(ds/dt)/ds)  half-life = log(2)/P(2) ~ 1sec
% P(2) - autoregulation   - d(ds/dt)/df)  2*pi*sqrt(1/P(3)) ~ 10 sec
% P(3) - transit time               (t0)  ~ 1 sec
% P(4) - exponent for Fout(v)    (alpha)  c.f. Grubb's exponent (~ 0.38)
% P(5) - resting oxygen extraction  (E0)  ~ range 20 - 50%

[qE,qC] = spm_hdm_priors(0);


% Reduce rank of prior covariances for computational expediancy
%---------------------------------------------------------------------------

% model specification (single node)
%---------------------------------------------------------------------------
M.fx  = 'spm_fx_dcm';
M.lx  = 'spm_lx_dcm';
M.x   = sparse(5,1);
M.m   = 1;
M.n   = 5;
M.l   = 1;
M.N   = 32;
M.dt  = 16/M.N;

% papermeters (C = 1)
%---------------------------------------------------------------------------
P     = [1; -1; 0; 1; qE];
p     = length(P);

% compute partial derivatives w.r.t. hemodynamic parameters [J] dy(t)/dp
%---------------------------------------------------------------------------
dp    = 1e-6;
M.pE  = P;
[k J] = spm_nlsi(M);
for i = 1:5
	M.pE    = P + sparse(4 + i,1,dp,p,1);
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


% combine connectivity and hemodynamic priors
%===========================================================================
qC    = kron(qC,eye(n,n));
qE    = kron(qE,ones(n,1));
pE    = [pE; qE];
pC    = blkdiag(pC,qC);


return


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
