
clear all
close all

disp('Model comparison using');
disp('Thermodynamic Sampling');

model_type='linear';
%model_type='logistic';

annealing=0;
model_switch=0;
ais=1;

switch model_type
    
    case 'linear',
        disp('Linear Regression');
        Nobs=20;
        sd=0.2;
        lambda=1/(sd^2);
        [M{1},U{1}] = mci_linear_struct (Nobs,lambda,'dct');
        w_true = spm_normrnd(M{1}.pE,M{1}.pC,1);
        Y = U{1}.X*w_true+sqrt(M{1}.Ce)*randn(Nobs,1);
        
        % Alternative model has fewer variables
        Nred=6; % Number of parameters in reduced model
        AM{1}=M{1};
        AU{1}.X=U{1}.X(:,1:Nred);
        AM{1}.pE=zeros(Nred,1);
        AM{1}.pC=M{1}.pC(1:Nred,1:Nred);
        
        % Indices of parameters to remove
        reduced_ind=[Nred+1:7];
        
    case 'logistic',
        disp('Logistic Regression');
        [M{1},U{1}] = mci_logistic_struct ('dct');
        P=[1 10 10]';
        [g,Y]=mci_logistic_gen(P,M{1},U{1});
        
        % Alternative model has no third variable
        reduced_ind=3;
        AM{1}=M{1};
        AU{1}.X=U{1}.X(:,1:2);
        AM{1}.pE=zeros(2,1);
        AM{1}.pC=100*eye(2);
        
end

% Set MCMC parameters
mcmc.nscale=500;
mcmc.ntune=500;
mcmc.nsamp=1000;
mcmc.J=16; % Number of temperatures/chains

mcmc.remove_burn_in=1;

% No sharing of samples between chains
mcmc.gprob=0;

if annealing
    disp('TI Annealing ...');
    tic;
    [P,logev1,D] = spm_mci_pop (mcmc,M,U,Y);
    toc
    tic;
    [P,logev2,D] = spm_mci_pop (mcmc,AM,AU,Y);
    toc
    disp(sprintf('TI Annealing Log Evidence Model 1 = %1.2f',logev1.ti));
    disp(sprintf('TI Annealing Log Evidence Model 2= %1.2f',logev2.ti));
end

if strcmp(model_type,'linear')
    [Ep,Cp,L1] = mci_linear_post (M{1},U{1},Y);
    [Ep,Cp,L2] = mci_linear_post (AM{1},AU{1},Y);
    disp(sprintf('Analytic log evidence Model 1 = %1.2f',L1));
    disp(sprintf('Analytic log evidence Model 2 = %1.2f',L2));
end

if model_switch
    % TI Model Switch
    M{2}=M{1};
    U{2}=U{1};
    for i=1:length(reduced_ind),
        k=reduced_ind(i);
        M{2}.pC(k,k)=1e-4;
    end
    mcmc.anneal='sigmoid';
    
    disp('TI Model Switching ...');
    tic;
    [Psd,logswitch,D] = spm_mci_pop (mcmc,M,U,Y);
    toc
end

if ais
    disp('Annealed Importance Sampling ...');
    mcmc.inference='ais';
    mcmc.anneal='geometric';
    mcmc.prop='lmc';
    mcmc.nprop=1;
    mcmc.J=512;
    mcmc.maxits=32;
    
    tic;
    post = spm_mci_ais (mcmc,M{1},U{1},Y);
    ais_logev1=post.logev;
    toc
    tic;
    post = spm_mci_ais (mcmc,AM{1},AU{1},Y);
    ais_logev2=post.logev;
    toc
    disp(sprintf('AIS Log Evidence Model 1 = %1.2f',ais_logev1));
    disp(sprintf('AIS Log Evidence Model 2= %1.2f',ais_logev2));
end

disp(' ');
if annealing
    disp(sprintf('TI Annealing Log Bayes Factor = %1.2f',logev2.ti-logev1.ti));
end
if strcmp(model_type,'linear')
    disp(sprintf('Analytic Log Bayes factor = %1.2f',L2-L1));
end
if model_switch
    disp(sprintf('TI Switching Log Bayes Factor = %1.2f',logswitch.ti));
end
if ais
    disp(sprintf('AIS Log Bayes Factor = %1.2f',ais_logev2-ais_logev1));
end