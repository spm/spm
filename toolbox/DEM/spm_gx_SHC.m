function [g] = spm_gx_SHC(x,v,P)
% maps to state of a SCH to a 2-D position in the world
% FORMAT [g] = spm_gx_SHC(x,v,P)
%
% x    - vector of hidden sates
% P.g  - state-space location associated with each hidden states
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_gx_SHC.m 3113 2009-05-11 15:25:13Z karl $
 
% expected position (coordinates in P.g)
%--------------------------------------------------------------------------
px     = exp(x(:));
px     = px/sum(px);
g      = P.g*px;
