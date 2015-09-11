function [M,U] = mci_rphase_init (d,conn)
% Initialise weakly coupled oscillator model - reduced connectivity
% FORMAT [M,U] = mci_rphase_init (d,conn)
%
% d     number of oscillators
%
% M     model structure
% U     inputs
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: mci_rphase_struct.m 6548 2015-09-11 12:39:47Z will $

M.N=100;
M.T=1;
dt=0.01;
M.t=[1:M.N]'*dt;
M.f='mci_rphase_fx';  
M.g='mci_phase_gx';
M.dfdp='mci_rphase_dfdp';
M.dfdx='mci_rphase_dfdx';

% Initial states
M.x0=rand(d,1);
M.x=zeros(d,1);
M.m=0;
M.n=d;
M.l=d;

nodes=M.n;

Nc=length(conn);
Np=2*Nc;

M.freq=6;
M.pE.aconn=zeros(Nc,1);
M.pE.bconn=zeros(Nc,1);
M.vpE=spm_vec(M.pE);
M.pC=eye(Np);

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

M.conn=conn;
