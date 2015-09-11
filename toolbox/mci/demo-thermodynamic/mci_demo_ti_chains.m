
clear all
close all

disp('Compute Model Evidence for Logistic Regression models using');
disp('Thermodynamic Integration');
disp('Explore effect of number of chains');

% Logistic model
[M{1},U{1}] = mci_logistic_struct ('dct');
P=[1 10 10]';
[g,Y]=mci_logistic_gen(P,M{1},U{1});
figure;
plot(g,'k','LineWidth',2);
set(gca,'FontSize',16);
grid on
xlabel('Time,t');
ylabel('p(Y=1|t)');

% MCMC parameters
mcmc.J=64; % Number of temperatures/chains
mcmc.ntune=500;
mcmc.nsamp=250;
setinit=1;
if setinit
    % Start sampling from this point
    for j=1:mcmc.J,
        mcmc.init{j}=P;
    end
    mcmc.nscale=0;
else
    mcmc.nscale=500;
end
mcmc.remove_burn_in=1;

% No sharing of samples between chains
mcmc.gprob=0;

j=[4,8,16,32,64];
for i=1:length(j),
    mcmc.J=j(i);
    tic;
    [P,logev,D] = spm_mci_pop (mcmc,M,U,Y);
    els(i)=toc;
    L(i)=logev.ti;
    disp(sprintf('J=%d: TI Log Evidence = %1.2f',j(i),logev.ti));
end

figure
plot(j,L,'k','LineWidth',2);
grid on
set(gca,'FontSize',16);
xlabel('Number of Chains');
ylabel('Log Evidence');
title('Thermodynamic Integration');

% Annealed Importance Sampling
mcmc.inference='ais';
mcmc.anneal='geometric';
mcmc.prop='lmc';
mcmc.nprop=1;
mcmc.maxits=64;


j=[4,8,16,32,64,128,256,512,1024];
for i=1:length(j),
    mcmc.J=j(i);
    tic;
    post = spm_mci_ais (mcmc,M{1},U{1},Y);
    els(i)=toc;
    L(i)=post.logev;
    disp(sprintf('J=%d: AIS Log Evidence = %1.2f',j(i),L(i)));
end

figure
semilogx(j,L,'k','LineWidth',2);
grid on
set(gca,'FontSize',16);
xlabel('Number of Chains');
ylabel('Log Evidence');
title('Annealed Importance Sampling');
