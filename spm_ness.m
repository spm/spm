function [H,R] = spm_ness(J,G)
% Evaluation of hessian and solenoidal operators at NESS
% FORMAT [H,R] = spm_ness(J,G)
% J  - Jacobian (dfdx)
% G  - diffusion tensor (amplitude of random fluctuations)
%
% H  - Hessian matrix (i.e., precision of a Gaussian density)
% R  - Skew symmetric solenoidal operator
%__________________________________________________________________________
% This routine evaluates the hessian (i.e., precision) of a nonequilibrium
% steady-state density (using a local linear approximation, under Gaussian
% assumptions). This is evaluated  under linear constraints on the
% solenoidal term of a Helmholtz decomposition. In short, given the flow
% (encoded by the systems Jacobian) and amplitude of random fluctuations,
% one can evaluate the steady-state density under nonequilibrium dynamics
% implied by solenoidal flow.
%
% flow constraints (Jacobian J)(R = -Q')
%--------------------------------------------------------------------------
% where flow   f = (R + G)*d log(p(x))/dx and
% log(p(x))      = -(1/2)*x'*H*x =>
% d log(p(x))/dx = -H*x =>
% df/dx = J      = -(R + G)*H =>
% H              = -(R + G)\J =>
% J*R + R*J'     = J*G - G*J'
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ness.m 7799 2020-03-12 17:23:14Z karl $


% solve for solenoidal (R) operator J*R + R*J' = J*G - G*J'
%==========================================================================
n = size(J,1);
I = eye(n,n);
X = kron(I,J) + kron(conj(J),I);
Y = spm_vec(J*G - G*J');
R = reshape(X\Y,n,n);

% precision (inverse covariance) of steady-state density
%--------------------------------------------------------------------------
H = -(R + G)\J;

return

% NOTES: stochastic coupling
%==========================================================================
clear HH RR
n     = 8;
for i = 1:64
    
    J = [1 1 1 0;
        1 1 1 0;
        0 1 1 1;
        0 1 1 1];
    J = kron(J,ones(n,n));
    J = J.*randn(size(J));
    J = J - eye(size(J))*4;
    G = diag(1 + rand(size(J,1),1));
    
    [H,R]     = spm_ness(J,G);
    
    RR(:,:,i) = R;
    HH(:,:,i) = H;
    JJ(:,:,i) = J;
    GG(:,:,i) = G;
end
subplot(2,2,1),imagesc(std(HH,0,3)), axis image, title('Hessian',   'Fontsize',16)
subplot(2,2,2),imagesc(std(RR,0,3)), axis image, title('Solenoidal','Fontsize',16)
subplot(2,2,3),imagesc(std(JJ,0,3)), axis image, title('Jacobian',  'Fontsize',16)
subplot(2,2,4),imagesc(std(GG,0,3)), axis image, title('Diffusion', 'Fontsize',16)


% complex analysis
%==========================================================================
[eJ,vJ] = eig(J,'nobalance','vector');
J       = pinv(eJ)*J*eJ;
G       = pinv(eJ)*G*eJ;

% complex form
%--------------------------------------------------------------------------
R = 1j*real(G)*imag(J)/real(J);
H = - real(J)/G;







