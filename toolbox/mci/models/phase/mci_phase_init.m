function [P,M,U,Y] = mci_phase_init (d)
% Initialise weakly coupled oscillator model
% FORMAT [P,M,U,Y] = mci_phase_init (d)
%
% d     number of oscillators
%
% P     parameters (drawn from prior)
% M     model structure
% U     inputs
% Y     data
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny and Biswa Sengupta
% $Id: mci_phase_init.m 6548 2015-09-11 12:39:47Z will $

M.N=200;
M.T=1;
dt=M.T/M.N;
M.t=[1:M.N]'*dt;
M.f='mci_phase_fx';  
M.g='mci_phase_gx';
M.dfdp='mci_phase_dfdp';
M.dfdx='mci_phase_dfdx';

% Initial states
M.x0=rand(d,1);
M.x=zeros(d,1);
M.m=0;
M.n=d;
M.l=d;

nodes=M.n;

params.cosCoeff = zeros(nodes);
params.sinCoeff = zeros(nodes);
params.intrinPhase = zeros(nodes,1);

pE=[zeros(2*nodes^2,1); 3*ones(nodes,1)];
Np=length(pE);
M.pE=spm_unvec(pE,params);
M.vpE=pE;
M.pC=0.1^2*eye(Np);
P=spm_normrnd(M.vpE,M.pC,1);

M.Np=Np;
M.ipC=inv(M.pC);

%M.int='euler';
%M.int='ode15';
M.int='sundials';
U=zeros(1,M.N);

% Likelihood function
obs_sd=0.1;
M.L='spm_mci_glike';
M.Ce=obs_sd^2*eye(d);

Y = spm_mci_fwd (P,M,U); 