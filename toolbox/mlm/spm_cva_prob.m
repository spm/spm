function [CVA] = spm_cva_prob (X1,X2,m)
% Probabilistic Canonical Variates Analysis
% FORMAT [CVA] = spm_cva_prob (X1,X2,m)
%
% X1           [d1 x N] matrix of dependent variables
% X2           [d2 x N] matrix of independent variables
% m            dimension of latent variable (min([d1,d2]) by default)
%
% Returns fields:
% 
% .U1,.U2      Canonical vectors
% .W1,.W2      Factor matrices
% .L           Log-Likelihood
% .bic         Bayesian Information Criterion
% .aic         Akaike's Information Criterion
%
% Fits probabilistic model
%
% x1 = W1 z + e1
% x2 = W2 z + e2
%
% This algorithm is described in:
%
% F. Bach and M. Jordan (2005) A probabilistic interpretation of canonical
% correlation analysis. Dept. Stats, Univ California, Berkeley CA. 
% Tech Rep 688.
%
%___________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_cva_prob.m 4687 2012-03-14 18:15:49Z will $


[d1,N1]=size(X1);
[d2,N2]=size(X2);
d=d1+d2;

if ~N1==N2
    disp('Error in spm_cva_prob: unequal number of samples');
    return
else
    N=N1;
end

if nargin < 3 | isempty(m)
    m=min([d1,d2]);
else
    if m > min([d1,d2]);
        disp('m too large');
        return
    end
end

X=[X1;X2];
Sigma=cov(X')+(10^-8*eye(d1+d2));
Sigma11=Sigma(1:d1,1:d1);

if m==0
    S=diag(diag(Sigma));
    iS=inv(S);
    CVA.L=-0.5*N*(d1+d2)*log(2*pi)-0.5*N*spm_logdet(S)-0.5*N*trace(iS*Sigma);
    CVA.bic=CVA.L;
    CVA.aic=CVA.L;
    CVA.W1=[];CVA.W2=[];CVA.U1=[];CVA.U2=[];
    return
end

Sigma12=Sigma(1:d1,d1+1:d1+d2);
Sigma22=Sigma(d1+1:d,d1+1:d);
Sigma21=Sigma12';

R1=inv(sqrtm(Sigma11));
R2=inv(sqrtm(Sigma22));
A=R1*Sigma12*R2;

[V1,P,V2]=svd(A,0);
sP=diag(P);
rP=sP/max(sP);
U1=R1*V1(:,1:m);
U2=R2*max(sP)*V2(:,1:m);

M1=diag(sqrt(rP(1:m)));
M2=M1;

W1=Sigma11*U1*M1;
W2=Sigma22*U2*M2;
Psi1=diag(diag(Sigma11-W1*W1'));
Psi2=diag(diag(Sigma22-W2*W2'));

W=[W1;W2];
S=W*W'+blkdiag(Psi1,Psi2);

CVA.W1=W1;
CVA.W2=W2;
CVA.U1=U1;
CVA.U2=U2;
 
iS=inv(S);
CVA.L=-0.5*N*(d1+d2)*log(2*pi)-0.5*N*spm_logdet(S)-0.5*N*trace(iS*Sigma);
k=2*m*(d1+d2); % W1, W2 plus diag Psi's 

CVA.bic=CVA.L-0.5*k*log(N);
CVA.aic=CVA.L-k;



