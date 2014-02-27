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
% $Id: spm_lx_ind.m 5900 2014-02-27 21:54:51Z karl $

% get lead field
%==========================================================================
G  = speye(M.l,M.l);


