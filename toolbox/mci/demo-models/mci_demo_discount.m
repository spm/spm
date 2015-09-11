
clear all
close all

disp('Temporal discounting model');

% Parameters for first level
Nobs=200;

[M,U] = mci_discount_struct (Nobs);

% Generate Data
%w_true = spm_normrnd(M.pE,M.pC,1);
w_true = [0.25, 1.5]';
w_true = log(w_true);
[g,Y] = mci_discount_gen (w_true,M,U);

[a,v1,v2,k] = mci_discount_act (w_true,M,U);

t=sort(U.t1);
r=mean(U.r1);
figure;plot(t,r./(1+k*t),'k.');
set(gca,'FontSize',16);
xlabel('Delay, Weeks');

t=sort(U.t2);
r=mean(U.r2);
hold on
plot(t,r./(1+k*t),'r.');
ylabel('Mean Reward, £');
grid on

disp(sprintf('Average prob of choosing first option = %1.2f',mean(g)));

%mcmc.inference='amc';
%mcmc.inference='vl';
mcmc.inference='langevin';

mcmc.verbose=0;
mcmc.maxits=1024;

post = spm_mci_post (mcmc,M,U,Y,w_true);
disp('Posterior Mean from MCI');
disp(post.Ep)

plot_post=1;
if plot_post
    % Plot posterior surface
    S.Nbins=50;
    r=2;
    S.pxy(1,:)=linspace(-r*w_true(1),r*w_true(1),S.Nbins);
    S.pxy(2,:)=linspace(-0.5*r*w_true(2),1.5*r*w_true(2),S.Nbins);
    S.param{1}='P(1)';
    S.param{2}='P(2)';
    S.name={'log k','log \beta'};
    [L,S] = mci_plot_surface (w_true,M,U,Y,S,'post');
    hold on
    ms=10;
    plot(M.pE(1),M.pE(2),'wo','MarkerSize',ms);
    plot(post.Ep(1),post.Ep(2),'wx','MarkerSize',ms);
    plot(w_true(1),w_true(2),'w+','MarkerSize',ms);
end

stats = spm_mci_mvnpost (post,'ESS')
stats = spm_mci_mvnpost (post,'thinning')

