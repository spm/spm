function [dLdp,iCpY] = spm_mci_grad_curve (assign,w,v,M,U,Y,fxtype)
% Compute gradient and curvature for MFX model
% FORMAT [dLdp,iCpY] = spm_mci_grad_curve (assign,w,v,M,U,Y,fxtype)
%
% assign    fields specify which are random/fixed effects
% w         random effects vector
% v         fixed effects vector
% M,U,Y     structure,inputs,data
% fxtype    'random' or 'fixed'
%
% dLdp      gradient
% iCpY      curvature
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_grad_curve.m 6275 2014-12-01 08:41:18Z will $

% Extract init and flow params from rfx or ffx vectors
[p_init,p_flow] = spm_mci_init_flow (assign,w,v,M);

% Gradient and curvature of log likelihood
if strcmp(assign.init_par,fxtype)
    [dLdp_init,iCpY_init] = spm_mci_grad_ind (p_flow,p_init,M,U,Y,'init');
else
    dLdp_init=[];iCpY_init=[];
end
if strcmp(assign.flow_par,fxtype)
    [dLdp_flow,iCpY_flow] = spm_mci_grad_ind (p_flow,p_init,M,U,Y,'flow');
else
    dLdp_flow=[];iCpY_flow=[];
end
dLdp=[dLdp_init, dLdp_flow];
iCpY=blkdiag(iCpY_init,iCpY_flow);

        