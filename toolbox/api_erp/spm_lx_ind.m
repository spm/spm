function [G] = spm_lx_ind(P,M)
% observer matrix for a DCM of induced responses: y = G*x
% FORMAT [G] = spm_lx_ind(P,M)
% x    - state vector - running over sources and then frequencies
%        NB: x(1) is time
% G    - y = G*x
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% %W% Karl Friston %E%

% parameterised lead field ECD
%--------------------------------------------------------------------------
L          = spm_erp_L(P,M);
Nf         = (M.n - 1)/M.r;
G          = sparse(M.l,M.n);
G(:,2:M.n) = kron(speye(Nf,Nf),L.^2);
