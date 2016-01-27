function [M,U,Y] = mci_logistic_struct (log_data,T)
% Set up data structures for logistic model
% FORMAT [M,U,Y] = mci_logistic_struct (log_data,T)
%
% log_data  'pima','ripley' or 'dct'
% T         for 'dct' we can specify number of samples
%
% M         model structure
% U         U.X is the design matrix (independent variables)
% Y         dependent variable
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_logistic_struct.m 6697 2016-01-27 14:57:28Z spm $

try, T=T; catch, T=100; end

switch log_data,
    % Data is from Brian Ripley's archive
    % http://www.stats.ox.ac.uk/pub/PRNN/
    % See also mci/data/README
    
    case 'pima',
        var_name={'npreg','glu','bp','skin','bmi','ped','age'};
        load pima_tr.txt
        load pima_te.txt
        U.X=[pima_tr(:,1:7);pima_te(:,1:7)];
        Y=[pima_tr(:,8);pima_te(:,8)];
        T=size(U.X,1);
        U.X=[U.X,ones(T,1)];
        prior_var=1;

    case 'ripley',
        var_name={'x1','x2'};
        load synth_tr.txt
        U.X=synth_tr(:,1:2);
        Y=synth_tr(:,3);        
        T=size(U.X,1);
        U.X=[U.X,ones(T,1)];
        prior_var=1;
        
    case 'dct',
        var_name={'x1','x2','x3'};
        U.X=spm_dctmtx(T,3);
        U.names=var_name;
        prior_var=100;
        Y=[];
        
    case 'cos',
        var_name={'x1','x2','x3'};
        U.names=var_name;
        Tmax=1000;
        t=[1:Tmax]/Tmax;
        U.X(:,1)=ones(Tmax,1);
        U.X(:,2)=cos(2*pi*t);
        U.X(:,3)=cos(4*pi*t);
        ind=randperm(Tmax);
        s=sort(ind(1:T));
        U.X=U.X(s,:);
        U.t=t(s);
        prior_var=100;
        Y=[];
end


M.L='mci_logistic_like';
M.dL='mci_logistic_deriv';
M.IS='mci_logistic_gen';
M.l=1;
M.T=T;
M.N=T;
Np=size(U.X,2);
M.pE=zeros(Np,1);
M.pC=prior_var*eye(Np);

