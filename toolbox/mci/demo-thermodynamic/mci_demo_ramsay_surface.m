
clear all
close all

disp('Synthetic data from nonlinear ODE model');
disp('defined in Ramsay et al. (2007)');
disp('based on Van der Pol oscillator and which');
disp('reduces to Fitzhugh-Nagumo for certain parameters');

disp(' ');
disp('Create data for posterior surface plot');

sigma_e=0.01;
[M,U] = mci_ramsay_struct(sigma_e);

% Use higher integration tolerances than by default
M.reltol=1e-3;
M.abstol=1e-5;
    
rand_init=0;
if rand_init
    P=spm_normrnd(M.vpE,M.pC,1);
else
    P=[log(0.2) log(0.2)]';
end

Y = mci_ramsay_gen(P,M,U);
mci_plot_outputs(M,Y);

S.Nbins=50;
S.pxy(1,:)=linspace(-2,0,S.Nbins);
S.pxy(2,:)=linspace(-2,0,S.Nbins);
S.param{1}='P(1)';
S.param{2}='P(2)';
S.name={'log a','log b'};
%[L,S] = mci_plot_surface (P,M,U,Y,S,'prior');
%[L,S] = mci_plot_surface (P,M,U,Y,S,'like');
tic;
[L,S] = mci_plot_surface (P,M,U,Y,S,'post');
toc

save ramsay-surface L S