function [M,U] = mci_nmm_struct (back,sd,Np)
% Set up two region NMM 
% FORMAT [M,U] = mci_nmm_struct (back,sd,Np)
%
% back      1 to include backward connection (default)
% sd        Observation noise SD (default 0.01)
% Np        number of params (2,6 or 21)
%
% M         Model structure
% U         Inputs
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: mci_nmm_struct.m 6697 2016-01-27 14:57:28Z spm $

if nargin < 1 || isempty(back)
    back=1; 
end
try, lambda=1/(sd^2); catch, lambda=10000; end

M.Nr=2;

%M.N=1000;
M.N=250;
M.T=0.4;
M.x=zeros(M.Nr,9);
M.x0=zeros(M.Nr*9,1);

M.n=length(M.x0);   % n states
M.l=2; % l outputs

% Specify input
M.ons=64;
M.dur=16;
PP.R=[0 0];
dt=M.T/M.N;
M.t=[1:M.N]'*dt;

M.g='mci_nmm_r2_gx';

% Specify connectivity pattern
A{1}=[0 0; 1 0]; % Extrinsic forward
if back
    A{2}=[0 1; 0 0]; % Extrinsic backward
else
    A{2}=[0 0; 0 0]; % Extrinsic backward
end
A{3}=[0 0; 0 0]; % Extrinsic lateral
B=[];

% Priors
% spm_erp_priors.m has prior SD of 0.25 on connex
prior_sd=0.4;
switch Np
    case 2,
        M.pE=zeros(2,1);
        M.pC=prior_sd^2*eye(2);
        M.m=2; % m inputs
        u(1,:)=spm_erp_u(M.t,PP,M);
        u(2,:)=zeros(1,M.N);
        U=full(u);
        M.f='mci_nmm_r2p2_fx';
        C=1;
        [M.can_P,pC] = spm_erp_priors(A,B,C);
        M.can_P.C=0;
        
    case 6,
        M.pE=zeros(6,1);
        M.pC=prior_sd^2*eye(6);
        M.m=2; % m inputs
        u(1,:)=spm_erp_u(M.t,PP,M);
        u(2,:)=zeros(1,M.N);
        U=full(u);
        M.f='mci_nmm_r2p6_fx';
        C=1;
        [M.can_P,pC] = spm_erp_priors(A,B,C);
        M.can_P.C=0;
    case 21,
        M.pC.A{1}(2,1)=prior_sd^2; % Forward connection
        M.pC.A{2}(1,2)=prior_sd^2; % Backward connection
        M.m=1; % m inputs
        u(1,:)=spm_erp_u(M.t,PP,M);
        U=full(u);
        M.f='mci_nmm_fx_delay';
        C=[1 0]';
        [M.pE,M.pC] = spm_erp_priors(A,B,C);
        M.vpE=spm_vec(M.pE);
    otherwise
        disp('Unknown Np value in nmm_struct.m');
end

%M.int='euler';
%M.int='ode15';
M.int='sundials';

%M.ipC=inv(M.pC);
M.Np=length(M.pE);
M.vpE=spm_vec(M.pE);

% Likelihood function
obs_sd=sqrt(1/lambda);
M.L='spm_mci_glike';
M.Ce=obs_sd^2*eye(M.Nr);

M.Npflow=length(M.pE);
M.Npout=0;

