
clear all
close all

disp('Synthetic data from d-region weakly coupled oscillator');
disp('Model with reduced (not full) connectivity');

% Dimension of state space
d=3;

conn(1).regions=[1,2];
conn(1).name='vis2mtl';
conn(2).regions=[3,2];
conn(2).name='vis2ifg';

[M,U] = mci_rphase_struct(d,conn);

M.freq=6;
M.x0=[-pi/2,0,pi/2]';

rand_init=0;
if rand_init
    P=spm_normrnd(M.vpE,M.pC,1);
else
    %pars.aconn=[2.75,-2.65]';
    %pars.bconn=[2.59,1.85]';
    pars.aconn=[-0.6 -0.6];
    pars.bconn=[0,0]';
    P=spm_vec(pars);
end

Y = mci_rphase_gen(P,M,U);

mci_plot_outputs(M,cos(Y));

%mcmc.inference='vl';
mcmc.inference='langevin';

mcmc.maxits=1024;
mcmc.verbose=0;

post = spm_mci_post (mcmc,M,U,Y,P);

spm_mci_diag(post);

stats = spm_mci_mvnpost (post,'ESS')
stats = spm_mci_mvnpost (post,'thinning')
stats = spm_mci_mvnpost (post,'MAR')


for j=1:length(P),
    spm_mci_quantiles (post,j,0);
end