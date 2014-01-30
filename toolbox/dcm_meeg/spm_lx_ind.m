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
% $Id: spm_lx_ind.m 5863 2014-01-30 20:58:36Z karl $

% get lead field
%==========================================================================
Nr         = length(P.L);           % number of sources
Nm         = M.n/Nr;                % number of frequency modes modelled
Nf         = M.l/Nr;                % number of frequency modes explained
G          = kron(speye(Nf,Nm),128*diag(exp(P.L)));


