
clear all
close all

disp('Linear Regression example');

% Design type
%des='linear-offset';
%des='quadratic-offset';
des='dct';
%des='quadratic';

% Parameters for first level
Nobs=20;
sd=0.2;
lambda=1/(sd^2);

[M,U,Xfull] = mci_linear_struct (Nobs,lambda,des);

% Generate Data
w_true = spm_normrnd(M.pE,M.pC,1);
Y = U.X*w_true+sqrt(M.Ce)*randn(Nobs,1);

%mcmc.inference='amc';
%mcmc.inference='vl';
mcmc.inference='langevin';

mcmc.verbose=0;
mcmc.maxits=1024;

post = spm_mci_post (mcmc,M,U,Y,w_true);
disp('Posterior Mean from MCI');
disp(post.Ep)

disp('Analytic Posterior Mean:');
[Ep,Cp]=mci_linear_post(M,U,Y);
disp(Ep)

dist{1}=post;  
dist{1}.color='k';
dist{1}.names=U.names;

dist{2}.Ep=Ep;
dist{2}.Cp=Cp;
dist{2}.names=dist{1}.names;
dist{2}.type='gaussian';
dist{2}.color='r';
figure
Np=length(M.pE);
for j=1:Np,
    subplot(1,Np,j);
    mci_plot_dist(dist{1},j);
    hold on
    mci_plot_dist(dist{2},j);
    grid on
end

ind=[1:size(post.Ce,3)];
mci_plot_noiseSD (post.Ce,ind);

% Plot model fits
names={'True','MCI','Analytic'};
y{1}=Xfull*w_true;
y{2}=Xfull*post.Ep;
y{3}=Xfull*Ep;
figure
mci_linear_plot_fit (M,Y,y,names);

spm_mci_diag(post);

stats = spm_mci_mvnpost (post,'ESS')
stats = spm_mci_mvnpost (post,'thinning')

[pstat,mu,nse,batch]=spm_mci_stat(post)

figure
plot(sqrt(squeeze(post.Ce(1,1,:))));
grid on
title('Noise SD');