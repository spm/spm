
clear all
close all

Nobs=40;

[M,U] = mci_approach_struct (Nobs);

rand_params=1;
if rand_params
    P = spm_normrnd(M.pE,M.pC,1);
else
    V = 30;
    tau = 8;
    P = [log(V),log(tau)]';
end

yhat = mci_approach_gen (P,M,U);
Y = yhat + sqrt(M.Ce)*randn(Nobs,1);

% Check analytic gradients using finite differences
%dLdp = approach_deriv (P,M,U,Y);
%dLdp = spm_diff(M.L,P,M,U,Y,1);

% figure
% plot(U.X,Y);
% hold on
% plot(U.X,yhat,'r');

%mcmc.inference='vl';
mcmc.inference='langevin';

mcmc.verbose=0;
mcmc.maxits=1024;

post = spm_mci_post (mcmc,M,U,Y,P);

disp('Prior Mean (circle):');
disp(M.pE);

disp('Posterior Mean from MCI (cross)');
disp(post.Ep)

disp('True (plus)');
disp(P)

% Plot posterior surface
S.Nbins=100;
S.pxy(1,:)=linspace(2,4,S.Nbins);
S.pxy(2,:)=linspace(0.5,3.5,S.Nbins);
S.param{1}='P(1)';
S.param{2}='P(2)';
S.name={'log V_a','log \tau'};
%[L,S] = mci_plot_surface (P,M,U,Y,S,'prior');
%[L,S] = mci_plot_surface (P,M,U,Y,S,'like');
[L,S] = mci_plot_surface (P,M,U,Y,S,'post');
hold on
ms=10;
plot(M.pE(1),M.pE(2),'wo','MarkerSize',ms);
plot(post.Ep(1),post.Ep(2),'wx','MarkerSize',ms);
plot(P(1),P(2),'w+','MarkerSize',ms);
j=post.ind;
plot(post.P(1,j),post.P(2,j),'w.');

diag.essplot=1
diag.ind=j;
mess=spm_mci_diag(post,diag);

stats = spm_mci_mvnpost (post,'ESS')
stats = spm_mci_mvnpost (post,'thinning')

