
clear all
close all

Nobs=40;

[M,U] = mci_pb_struct (Nobs);
M.Ce=4;

P = spm_normrnd(M.pE,M.pC,1);

yhat = mci_pb_gen (P,M,U);
Y = yhat + sqrt(M.Ce)*randn(Nobs,1);

%mcmc.inference='vl';
mcmc.inference='langevin';

mcmc.verbose=0;
mcmc.maxits=2048;

post = spm_mci_post (mcmc,M,U,Y,P);

disp('Prior Mean, Posterior Mean, True:');
disp([M.pE,post.Ep,P]);

y_true = mci_pb_gen (P,M,U);
y_prior = mci_pb_gen (M.pE,M,U);
y_post = mci_pb_gen (post.Ep,M,U);

figure
plot(U.X,y_true,'b');
hold on
plot(U.X,y_post,'r');
plot(U.X,y_prior,'k');
grid on
plot(U.X,Y,'.');

j=post.ind;
diag.essplot=1;
diag.ind=j;
mess=spm_mci_diag(post,diag);

disp('p-value for test of multivariate normality:');
p = spm_mci_mvntest(post.P(:,j)',mess)

