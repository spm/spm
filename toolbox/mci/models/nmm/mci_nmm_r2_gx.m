function [y,L] = mci_nmm_r2_gx (x,u,P,M)
% Observation function for 2-region NMM
% FORMAT [y,L] = mci_nmm_r2_gx (x,u,P,M)
%
% P         Parameters
% M         Model structure
% U         Inputs
%
% y         Output
% L         Lead field (dy/dx)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: mci_nmm_r2_gx.m 6548 2015-09-11 12:39:47Z will $

x=spm_vec(x);

L=zeros(2,18);
L(1,17)=1;
L(2,18)=1;
y=L*x;



