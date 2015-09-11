function [pinit,pflow,names,M,U,Y] = mci_lds_group_data (lds)
% Generate LDS data for a group of subjects
% FORMAT [pinit,pflow,names,M,U,Y] = mci_lds_group_data (lds)
%
% lds        Data structure with fields:
%
% .R         R.pE, R.pC prior over initial conds
% .sd        Standard deviation of observation noise
% .Nsub      Number of subjects
% .Nobs      Number of observations per subject
% .model     'lds_real','forward',etc.
% .flow_par  'fixed' or 'random'
% .init_par  'fixed' or 'random'
% 
% pinit      Initial params
% pflow      Flow params
% names      names of parameters
% M          Cell of models
% U          Cell of inputs
% Y          Cell of data
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_lds_group_data.m 6548 2015-09-11 12:39:47Z will $

R=lds.R; sd=lds.sd; Nsub=lds.Nsub; 
Nobs=lds.Nobs; model=lds.model;
d=length(R.pE);

% Set Flow Parameters
if strcmp(model,'lds_real')
    r=linspace(0.3,0.9,d);
    T=100;
    q=(1/T)*log(r);
    pflow=log(-q);
else
    
    Mtmp.name=model;
    Mtmp.sd=sd;
    Mtmp.d=d;
    Mtmp.drop=0.5;
    %Mtmp.t=[1:400]'/4; 
    Mtmp.t=[1:25]'/0.25; 
    Mtmp.int='sundials';
    
    Mtmp=mci_lds_struct(Mtmp);
    if strcmp(lds.flow_par,'fixed')
%         PP.self=Mtmp.sd_self*randn(d,1);
%         PP.between=Mtmp.sd_between*randn(Mtmp.Nb,1);
%         pflow=spm_vec(PP);
        
        Pt = [-0.04,-0.01,-0.005,-0.01,0.01,0.01,0.01]';
        pflow = mci_lds_par2lat (Pt,Mtmp);
    else
        for n=1:Nsub,
            PP.self=Mtmp.sd_self*randn(d,1);
            PP.between=Mtmp.sd_between*randn(Mtmp.Nb,1);
            pflow(:,n)=spm_vec(PP);
        end
    end
end

if strcmp(lds.init_par,'fixed')
    R0 = spm_normrnd(R.pE,R.pC,1);
end
for n=1:Nsub,
    
    % Sample initial states from prior
    if strcmp(lds.init_par,'random')
        R0 = spm_normrnd(R.pE,R.pC,1);
    end
    pinit(:,n)=R0;
    switch model
        case 'lds_real'
            [M{n},U{n},y] = irlds_init (d,sd,R0,pflow);
        otherwise
            Mtmp.R=R0;
            [M{n},U{n}] = mci_lds_struct (Mtmp);
            
            if strcmp(lds.flow_par,'fixed')
                y = mci_lds_gen (M{n},U{n},pflow);
            else
                y = mci_lds_gen (M{n},U{n},pflow(:,n));
            end
    end
    
    if lds.Nobs==M{n}.N
        Y{n}.y=y;
        Y{n}.ind=1:M{n}.N;
    else
        % Thin observations to selected time points
        rind=randperm(M{n}.N);
        ind=rind(1:Nobs);
        Y{n}.y=y(ind,:);
        Y{n}.ind=ind;
    end
    
end

if strcmp(model,'forward')
    
    for j=1:d,
        jn=int2str(j);
        names{j}=['a_{',jn,jn,'}'];
    end
    for i=1:d-1,
        jn=int2str(i);
        j1n=int2str(i+1);
        names{i+j}=['a_{',j1n,jn,'}'];
    end
    
end
