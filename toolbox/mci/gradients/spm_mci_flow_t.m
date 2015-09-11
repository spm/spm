function [dxdt] = spm_mci_flow_t (t,x,U,P,M)
% Evaluate flow at time t
% FORMAT [dxdt] = spm_mci_flow_t (t,x,U,P,M)
%
% t     time
% x     state
% U     inputs
% P     parameters
% M     model
%
% dxdt  flow, dx/dt
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: spm_mci_flow_t.m 6548 2015-09-11 12:39:47Z will $

% Find nearest time point for which we have pre-computed input
% (We could also compute input on the fly)
if isempty(U)
    ut=[];
else
    [tmp,ind]=min(abs(t-M.t));
    ut=U(:,ind);
end

dxdt=feval(M.f,x,ut,P,M);
