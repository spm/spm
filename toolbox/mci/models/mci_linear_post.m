function [Ep,Cp] = mci_linear_post (M,U,Y)
% Analytic posterior for linear regression
% FORMAT [Ep,Cp] = mci_linear_post (M,U,Y)
% 
% M     Model Structure
% U     Inputs
% Y     Data
%
% Ep    Posterior mean
% Cp    Posterior covariance
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_linear_post.m 6275 2014-12-01 08:41:18Z will $

ipC=inv(M.pC);
iCe=inv(M.Ce);
X=U.X;
iCp=X'*iCe*X+ipC;
Cp=inv(iCp);
Ep=Cp*(X'*iCe*Y+ipC*M.pE);