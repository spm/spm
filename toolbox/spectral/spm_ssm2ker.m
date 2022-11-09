function [ker,pst] = spm_ssm2ker(dfdx,dfdu,dgdx,pst)
% computes cross spectral density from state space representation
% FORMAT [ker,pst] = spm_ssm2ker(dfdx,dfdu,dgdx,pst)
%
% dfdx - Jacobian
% dfdu - input matrix  [default: 1]
% dgdx - output matrix [default: 1]
% pst  - time          [default: based on maximum eigenvalue]
%
% ker  - first-order (Volterra) kernels
% 
% NB: Please see notes at the end of this routine for a demonstration of
% the systems analyses using the suite of spm_???2??.m routines
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% preliminaries
%--------------------------------------------------------------------------
if nargin < 2,    dfdu = [];                end
if nargin < 3,    dgdx = [];                end
if isempty(dfdu), dfdu = speye(size(dfdx)); end
if isempty(dgdx), dgdx = speye(size(dfdx)); end

% Jacobian and eigenspectrum
%==========================================================================
try
   [v,s] = eig(full(dfdx),'nobalance');
catch
   v  = eye(size(dfdx));
   s  = NaN(size(dfdx));
end
s     = diag(s);

% condition unstable eigenmodes
%--------------------------------------------------------------------------
s = 1j*imag(s) + min(real(s),-1/64);

% time lags if unspecified
%--------------------------------------------------------------------------
if nargin < 4
    pst   = -8/max(real(s));
    pst   = pst*(1:64)'/64;
end

% kernel functions
%==========================================================================
ng    = size(dgdx,1);         % number of outputs
nu    = size(dfdu,2);         % number of inputs
nk    = size(v,2);            % number of modes
nt    = size(pst(:),1);       % number of time points
ker   = zeros(nt,ng,nu);      % kernels

% derivatives over modes
%--------------------------------------------------------------------------
dgdv  = dgdx*v;
dvdu  = pinv(v)*dfdu;
for j = 1:nu
    for i = 1:ng
        for k = 1:nk
            
            % kernels
            %--------------------------------------------------------------
            Kk         = exp(s(k)*pst(:));
            ker(:,i,j) = ker(:,i,j) + real(dgdv(i,k)*dvdu(k,j)*Kk);
            
        end
    end
end

return


% NOTES: for linear systems analysis
%==========================================================================

% create nearly critical Jacobian
%--------------------------------------------------------------------------
s     = 1;
while s > -1/4 || s < -1/2
    J = randn(8,8) - eye(8,8)*2;
    e = eig(J);
    s = max(real(e));
end

%  creat a dynamic causal model and solve for a realised timeseries
%--------------------------------------------------------------------------
M.f   = @(x,u,P,M) P*x;
M.x   = randn(8,1)/32;
U     = sparse(1024,8);
y     = spm_int_sde(J,M,U);

% linearise and evaluate transfer functions
%--------------------------------------------------------------------------
dfdx       = spm_dcm2ssm(J,M);
[mtf,Hz]   = spm_ssm2mtf(dfdx);
[ker,pst]  = spm_ssm2ker(dfdx);
dt         = pst(2) - pst(1);
dw         =  Hz(2) -  Hz(1);

% alternative evaluations of cross covariance function
%--------------------------------------------------------------------------
[ccf1,pst1] = spm_mtf2ccf(mtf,Hz);
[ccf2,pst2] = spm_ker2ccf(ker,dt);
[ccf3,pst3] = spm_ssm2ccf(dfdx,[],[],Hz);

% predicted and observed correlations
%--------------------------------------------------------------------------
i           = 1;
j           = 4;
[xc,lags]   = xcov(y(:,i),y(:,j),numel(pst) - 1,'biased');

% plot results
%==========================================================================
subplot(3,2,1), plot(real(y(:,[i,j])))
title('Time series'), xlabel('time'), axis square

subplot(3,2,2), plot(pst,real(ker(:,i,j)),pst,imag(ker(:,i,j)),':')
title('Kernels'), xlabel('time'), axis square

subplot(3,2,3), plot(Hz,real(mtf(:,i,j)),Hz,imag(mtf(:,i,j)),':')
title('Modulation transfer functions'), xlabel('frequency'), axis square

subplot(3,2,4), plot(pst1,ccf1(:,i,j),pst2,ccf2(:,i,j),pst3,ccf3(:,i,j),lags*dt,xc)
title('Covariance functions'), xlabel('time'), axis square
legend({'mtf2ccf','ker2ccf','ssm2ccf','sample'}), legend('boxoff')

subplot(3,2,5), imagesc(spm_ccf2cor(ccf3))
title('Predicted correlations'), axis square

subplot(3,2,6), imagesc(corr(y))
title('Realised correlations'),  axis square

