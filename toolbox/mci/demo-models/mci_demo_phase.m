
clear all
close all

disp('Synthetic data from d-region weakly coupled oscillator');
disp('Fully connected network');

M.N=200;
M.T=1;
dt=M.T/M.N;
M.t=[1:M.N]'*dt;
M.f='mci_phase_fx';  
M.g='mci_phase_gx';

% Dimension of state space
d=3;
[P,M,U,Y] = mci_phase_init(d);

mci_plot_outputs(M,cos(2*pi*Y));

% Check to see if VL will run
%tmp=spm_int(M.pE,M,U);

%mcmc.inference='vl';
mcmc.inference='langevin';

mcmc.verbose=1;
mcmc.maxits=256;

post = spm_mci_post (mcmc,M,U,Y,P);