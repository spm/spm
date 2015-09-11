function [Ep,Cp,L] = mci_linear_post (M,U,Y)
% Analytic posterior for linear regression
% FORMAT [Ep,Cp,L] = mci_linear_post (M,U,Y)
% 
% M     Model Structure
% U     Inputs
% Y     Data
%
% Ep    Posterior mean
% Cp    Posterior covariance
% L     Log evidence
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_linear_post.m 6548 2015-09-11 12:39:47Z will $

ipC=inv(M.pC);
iCe=inv(M.Ce);
X=U.X;
iCp=X'*iCe*X+ipC;
Cp=inv(iCp);
Ep=Cp*(X'*iCe*Y+ipC*M.pE);

% Log evidence
yhat=X*Ep;
T=length(yhat);
if isstruct(Y)
    ey=Y.y-yhat;
else
    ey=Y-yhat;
end

L = -0.5*T*spm_logdet(M.Ce) - 0.5*T*log(2*pi);
L = L - 0.5*trace(ey'*iCe*ey);

ew = M.pE-Ep;
L = L - 0.5*ew'*ipC*ew - 0.5*spm_logdet(M.pC) + 0.5*spm_logdet(Cp);


% Model evidence
% e1=y-glm.yhat;
% e2=glm.mu-glm.Ep;
% accuracy=-0.5*N*log(glm.c1)-0.5*(e1'*e1)/glm.c1;
% complexity=0.5*p*log(glm.c2)-0.5*log(det(glm.Cp))+0.5*(e2'*e2)/glm.c2;
% const=-0.5*N*log(2*pi);
% glm.logev=accuracy-complexity+const;
% glm.like=accuracy+const;
% glm.acc=accuracy;
% glm.kl=complexity;