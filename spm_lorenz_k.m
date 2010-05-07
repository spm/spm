function [f] = spm_lorenz_k(x,v,P)
% Equations of motion for coupled Lorenz attractors
% FORMAT [f] = spm_lorenz_k(x,v,P)
% x - hidden states (3 x N)
% v - exogenous input
% P - parameters
%     P.t = N x 1
%     P.k = 1 x 1
%__________________________________________________________________________
 
% take expectations and unpack
%--------------------------------------------------------------------------
X  = mean(x,2);                            % order parameters (mean-field)
N  = size(x,2);                            % number of oscillators
 
% p(1): Prandtl number
% p(2): 8/3
% p(3): Rayleigh number
%--------------------------------------------------------------------------
p  = [10; -8/3; 32];
f  = x;
for i = 1:N
    f(:,i) = [-p(1) p(1) 0; p(3) -1 -x(1,i); x(2,i) 0 p(2)]*x(:,i);
    f(:,i) = (f(:,i) + P.k*X)*exp(P.t(i));
end
 
% vectorised flow
%--------------------------------------------------------------------------
f  = spm_vec(f);
