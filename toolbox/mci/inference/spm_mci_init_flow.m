function [x_init,x_flow] = spm_mci_init_flow (assign,w,v,M)
% Extract init, flow and out params from rfx and ffx vectors
% FORMAT [x_init,x_flow] = spm_mci_init_flow (assign,w,v,M)
%
% assign    fields specify which are random/fixed effects
% w         random effects vector
% v         fixed effects vector
% M         model structure
%
% x_init    init params
% x_flow    flow params (includes out params)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_init_flow.m 6548 2015-09-11 12:39:47Z will $

switch assign.init_par,
    case 'known',
        x_init=M.x0;
    case 'random'
        x_init=w(assign.w_init);
    otherwise
        % Assume fixed effect
        x_init=v(assign.v_init);
end

if strcmp(assign.flow_par,'random')
    x_flow=w(assign.w_flow);
else
    x_flow=v(assign.v_flow);
end

if strcmp(assign.out_par,'random')
    x_flow=[x_flow;w(assign.w_out)];
else
    x_flow=[x_flow;v(assign.v_out)];
end
