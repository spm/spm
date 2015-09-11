function [M,U,Xfull] = mci_linsqr_struct (Nobs,lambda,des)
% Set up data structures for linsqr model
% FORMAT [M,U,Xfull] = mci_linsqr_struct (Nobs,lambda,des)
%
% Nobs      number of data points
% lambda    noise precision
% des       type of design 
%
% M         model structure
% U         U.X is the design matrix
% Xfull     Design matrix for data points [1:T]
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_linsqr_struct.m 6548 2015-09-11 12:39:47Z will $

% Number of time points
T=100;

order=2;
for j=1:order,
    U.names{j}=sprintf('R%d',j);
end
U.X=spm_dctmtx(T,order);

Xfull=U.X;

t=[1:T]';
% Thin observations to selected time points
if Nobs < T
    rind=randperm(T);
    ind=rind(1:Nobs);
    ind=sort(ind);
    U.X=U.X(ind,:);
    t=t(ind);
    U.ind=ind;
end
    
sigma_e=sqrt(1/lambda);
M.Ce=1/lambda;
M.L='mci_linsqr_like';
M.dL='mci_linsqr_deriv';
M.IS='mci_linsqr_gen';

Np=size(U.X,2);
M.pE=zeros(Np,1);
%M.pC=0.5*eye(Np);
M.pC=100*eye(Np);
M.l=1;
M.t=t;
M.T=T;
M.N=T;


