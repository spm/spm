function [f,J,Q] = spm_fx_tfm(x,u,P,M)
% state equations - time-frequency model with state-dependent parameters
% FORMAT [f,J,D] = spm_fx_tfm(x,u,P,M)
% x      - hidden states
% u      - exogenous input
%
% f        - dx(t)/dt  = f(x(t))
% J        - df(t)/dx(t)
% D        - delay operator dx(t)/dt = f(x(t - d))
%                                    = D(d)*f(x(t))
%
% This routine is essentially a rapper for the equations of motion
% specified in M.h - it updates the input dependent parameters and then
% calls the appropriate equations of motion in the usual way.
%
%__________________________________________________________________________
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%___________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_fx_tfm.m 7679 2019-10-24 15:54:07Z spm $
 
% input and state-dependent parameters
%==========================================================================
if isfield(M,'u')
    
    % endogenous inputs
    %----------------------------------------------------------------------
    P   = rmfield(P,{'X','Y'});
    
else
    
    % exogenous inputs
    %----------------------------------------------------------------------
    dP  = P.X*u(:) + P.Y*x(:);
    P   = rmfield(P,{'X','Y'});
    P   = spm_unvec(spm_vec(P) + dP,P);
    
end


% Equations of motion - place in model and evaluate
%==========================================================================
M.f = M.h;

% and evaluate
%--------------------------------------------------------------------------
if nargout == 3
    [f,J,Q] = feval(M.f,x,u,P,M);
    
elseif nargout == 2
    [f,J]   = feval(M.f,x,u,P,M);
    
else
    f       = feval(M.f,x,u,P,M);
end
