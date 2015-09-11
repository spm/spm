
clear all
close all

disp('Thermodynamic Integration Example');

%model_type='linear';
%model_type='growth';
model_type='logistic-dct';

switch model_type
    
    case 'linear',
        Nobs=20;
        sd=0.2;
        lambda=1/(sd^2);
        [M{1},U{1}] = mci_linear_struct (Nobs,lambda,'dct');
        w_true = spm_normrnd(M{1}.pE,M{1}.pC,1);
        Y = U{1}.X*w_true+sqrt(M{1}.Ce)*randn(Nobs,1);
        
    case 'growth',
        Nobs=40;
        [M{1},U{1}] = mci_pb_struct (Nobs);
        M{1}.Ce=4;
        Ptrue = spm_normrnd(M{1}.pE,M{1}.pC,1);
        yhat = mci_pb_gen (Ptrue,M{1},U{1});
        Y = yhat + sqrt(M{1}.Ce)*randn(Nobs,1);
        
    case 'logistic-ripley',
        [M{1},U{1},Y] = mci_logistic_struct ('ripley');
        
    case 'logistic-pima';
        [M{1},U{1},Y] = mci_logistic_struct ('pima');
        
    case 'logistic-dct';
        [M{1},U{1},Y] = mci_logistic_struct ('dct');
        Ptrue=[1 10 10]';
        [g,Y]=mci_logistic_gen(Ptrue,M{1},U{1});
end

% Set MCMC parameters
mcmc.nscale=250;
mcmc.ntune=250;
mcmc.nsamp=512;
mcmc.J=64; % Number of temperatures/chains

mcmc.remove_burn_in=1;

% No sharing of samples between chains
mcmc.gprob=0;

tic;
[P,logev,D] = spm_mci_pop (mcmc,M,U,Y);
toc

disp(sprintf('TI Log Evidence = %1.2f',logev.ti));

% Annealed Importance Sampling
mcmc.inference='ais';
mcmc.anneal='geometric';
mcmc.prop='lmc';
mcmc.nprop=1;
mcmc.J=512;
mcmc.maxits=64;

tic;
post = spm_mci_ais (mcmc,M{1},U{1},Y);
toc
disp(sprintf('AIS Log Evidence = %1.2f',post.logev));

if strcmp(model_type,'linear')
    [Ep,Cp,L] = mci_linear_post (M{1},U{1},Y);
    disp(sprintf('Analytic log evidence = %1.2f',L));
end