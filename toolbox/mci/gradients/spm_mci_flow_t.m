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

% Will Penny and Biswa Sengupta
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Find nearest time point for which we have pre-computed input
% (We could also compute input on the fly)
if isempty(U)
    ut=[];
else
    [tmp,ind]=min(abs(t-M.t));
    ut=U(:,ind);
end

dxdt=feval(M.f,x,ut,P,M);
