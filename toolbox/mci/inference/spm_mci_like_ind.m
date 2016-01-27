function [L,e] = spm_mci_like_ind (P,R,M,U,Y)
% Compute likelihood wrt selected time points
% FORMAT [L,e] = spm_mci_like_ind (P,R,M,U,Y)
%
% P         Flow parameters
% R         Initial state parameters
% M         Model structure
% U         Inputs  [Nin x N]
% Y         data
%     
% L         Log likelihood
% e         Prediction errors
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: spm_mci_like_ind.m 6697 2016-01-27 14:57:28Z spm $

% Read data points and time indices
try, ind=Y.ind; catch, ind=1:M.N; end
Nt=length(ind);
y=Y.y;

M.x0=R; % Initial conditions
[G,sy,st] = spm_mci_fwd (P,M,U);

if st==-1, disp('Problem !'); return; end

% Prediction errors
g=G(ind,:);
e=Y.y-g;

% Log Likelihood
L = -0.5*trace(M.iCe*e'*e) + M.logdet_Ce - 0.5*Nt*log(2*pi);