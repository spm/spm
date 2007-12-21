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
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_lx_ind.m 1040 2007-12-21 20:28:30Z karl $

% get lead field
%==========================================================================
Nr         = length(P.L);                      % number of sources
Nf         = (M.n - 1)/Nr;                     % number of frequency modes
L          = kron(speye(Nf,Nf),128*diag(exp(P.L)));
G          = sparse(M.l,M.n);
G(:,2:M.n) = L;

