function [G] = spm_lx_ind(P,M)
% observer matrix for a DCM of induced responses: y = G*x
% FORMAT [G] = spm_lx_ind(P,M)
% x    - state vector - running over sources and then frequencies
%        
% G    - y = G*x
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_lx_ind.m 1174 2008-02-27 20:22:30Z karl $

% get lead field
%==========================================================================
Nr         = length(P.L);                      % number of sources
Nf         = M.n/Nr;                           % number of frequency modes
G          = kron(speye(Nf,Nf),128*diag(exp(P.L)));


