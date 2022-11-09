function [f,dfdx,dfdu,dfdxdu] = spm_fx_mfm_D(x,u,P,M)
% State equations for neural-mass and mean-field models (delay version)
% FORMAT [f,dfdx,dfdu,dfdxdu] = spm_fx_mfm(x,u,P,M)
%
% x         - states (means and covariances)
% dfdx,...  - derivatives with repect to x and u
%
% x{1}(i,j,k)   - k-th state of j-th population on i-th source
% x{2}(:,:,i,j) - covariance among k states
%
%   population: 1 - excitatory spiny stellate cells (input cells)
%               2 - inhibitory interneurons
%               3 - excitatory pyramidal cells      (output cells)
%
%        state: 1 V  - voltage
%               2 gE - conductance (excitatory)
%               3 gI - conductance (inhibitory)
%
%--------------------------------------------------------------------------
% refs:
%
% This routine is exactly the same as spm_fx_mfm but premultiplies the flow
% with the delay operator to return the flow on delayed states. This is
% necessary for accurate computation of the Jacobian under steady state
% assumptions
%
% Delays
%==========================================================================
% Delay differential equations can be integrated efficiently (but 
% approximately) by absorbing the delay operator into the Jacobian
%
%    f(d) = dx(t)/dt = f(x(t - d))
%                    = Q(d)f(x(t))
%
%    J(d)            = Q(d)df/dx
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging

 
% get delayed flow and Jacobian
%--------------------------------------------------------------------------
[f,J,Q] = spm_fx_mfm(x,u,P,M);
f       = Q*f;
J       = Q*J;
