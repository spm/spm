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


% get lead field
%==========================================================================
Nr         = length(P.L);                      % number of sources
Nf         = (M.n - 1)/Nr;                     % number of frequencies
L          = kron(speye(Nf,Nf),diag(P.L));
G          = sparse(M.l,M.n);
G(:,2:M.n) = L;

return




% NOTES from channel space: parameterised lead field ECD given positions
%--------------------------------------------------------------------------
M.dipfit.type  = 'ECD (EEG)';
pos  = M.dipfit.L.pos;
mom  = [1  0  0;
        0  1  0;
        0  0  1;
        1  1  0;
        1 -1  0;
        1  0  1;
        1  0 -1;
        0  1  1;
        0  1 -1;
        1  1 -1;
        1 -1  1;
       -1  1  1
        1  1  1]';
   
Ng     = size(mom,2);                  % number of moments per source
Nr     = size(pos,2);                  % number of sources
G.Lmom = kron(ones(1,Nr),mom);
G.Lpos = kron(pos,ones(1,Ng));
L      = spm_erp_L(G,M);

% Gain matrix (for power) - accounting for x(1) = time
%--------------------------------------------------------------------------
for i = 1:Nr
    g{i} = exp(spm_vec(exp(P.L(:,i))))/16;
end
L          = (L.^2)*spm_cat(diag(g));
G          = sparse(M.l,M.n);
Nf         = (M.n - 1)/Nr;
G(:,2:M.n) = kron(speye(Nf,Nf),L);

