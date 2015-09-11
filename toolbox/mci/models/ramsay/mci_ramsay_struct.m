function [M,U] = mci_ramsay_struct (sigma_e)
% Data structures for Ramsay model
% FORMAT [M,U] = mci_ramsay_struct (sigme_e)
%
% sigma_e       Noise SD
%
% M,U           model, input data structures
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_ramsay_struct.m 6548 2015-09-11 12:39:47Z will $

ms=20;
dt=0.1;

T=ms/dt;
U.dt=dt;
U.u=zeros(T,1);
U.tims=[1:T]'*dt;

M.t=U.tims;
M.T=T;
M.N=length(U.tims);

M.x0  = [-1 1]';
M.x = M.x0;
%M.ns  = T;

M.m=1; % number of inputs
M.n=length(M.x0);   % n states
M.l=2; % l outputs

M.f   = 'mci_ramsay_fx';
M.g   = 'mci_ramsay_gx';
M.L   = 'spm_mci_glike';

%M.int='euler';
%M.int='ode15';
M.int='sundials';

M.pE = [log(0.5),log(0.5)]';
M.vpE = spm_vec(M.pE);
M.pC = [1/8 0; 0 1/8];
M.ipC = inv(M.pC);
M.Np = length(M.pE);

M.Ce  = sigma_e^2*eye(2);
