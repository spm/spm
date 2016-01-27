function [M,U,Xfull] = mci_linear_struct (Nobs,lambda,des)
% Set up data structures for linear model
% FORMAT [M,U,Xfull] = mci_linear_struct (Nobs,lambda,des)
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
% $Id: mci_linear_struct.m 6697 2016-01-27 14:57:28Z spm $

try, design=des; catch, design='linear-offset'; end

% Number of time points
T=100;

switch design
    case 'dct',
        order=7;
        for j=1:order,
            U.names{j}=sprintf('R%d',j);
        end
        U.X=spm_dctmtx(T,order);
    case 'linear-offset',
        x1=ones(T,1);
        x2=[1:T]';
        x2=x2-mean(x2);
        U.X=[x1,x2];
        U.names={'Offset','Linear'};
    case 'quadratic-offset',
        x0=ones(T,1);
        x1=[1:T]';
        x1=x1-mean(x1);
        x2=x1.^2;
        U.X=[x0,x1,x2];
        U.names={'Offset','Linear','Quadratic'};
    case 'quadratic',
        x1=[1:T]';
        x1=x1-mean(x1);
        x2=x1.^2;
        U.X=[x1,x2];
        U.names={'Offset','Quadratic'};
    case 'unequal_var',
        x1=[1:T]';
        x1=x1-mean(x1);
        x2=x1.^2;
        x2=x2/10;
        U.X=[x2,x1];
        U.names={'Offset','Quadratic'};
    otherwise
        disp(sprintf('Unknown design: %s, in linear_struct.m',design));
        return
end

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
M.L='mci_linear_like';
M.dL='mci_linear_deriv';
M.IS='mci_linear_gen';

Np=size(U.X,2);
M.pE=zeros(Np,1);
%M.pC=0.5*eye(Np);
M.pC=10*eye(Np);
M.l=1;
M.t=t;
M.T=T;
M.N=T;

M.Npflow=Np;
M.Npout=0;


