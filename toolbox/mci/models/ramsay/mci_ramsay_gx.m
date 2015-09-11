function [y,L] = spm_ramsay_gx (x,u,P,M)
% Observation equation for Ramsay model
% FORMAT [y,L] = spm_ramsay_gx (x,u,P,M)
%
% x,u,P,M     state,input,params,model
%
% y           observations
% L           dy/dx
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: mci_ramsay_gx.m 6548 2015-09-11 12:39:47Z will $

nx=length(x);
L=eye(nx);
y=x;

