function [G] = spm_lx_ind(P,M)
% Observer matrix for a DCM of induced responses: y = G*x
% FORMAT [G] = spm_lx_ind(P,M)
% x    - state vector - running over sources and then frequencies
%        
% G    - y = G*x
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging

% get lead field
%==========================================================================
G  = speye(M.l,M.l);
