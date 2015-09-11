
clear all
close all

disp('Linear regression with squared parameters');

% Design type
%des='linear-offset';
%des='quadratic-offset';
des='dct';
%des='quadratic';

% Parameters for first level
Nobs=100;
sd=0.1;
lambda=1/(sd^2);

[M,U,Xfull] = mci_linsqr_struct (Nobs,lambda,des);

% Generate Data
%w_true = spm_normrnd(M.pE,M.pC,1);
w_true = [2, 2]';
Y_true = mci_linsqr_gen (w_true,M,U);
Y = Y_true + sqrt(M.Ce)*randn(Nobs,1);

% figure;
% lw=2;
% plot(M.t,Y_true,'k','LineWidth',lw);
% hold on
% plot(M.t,Y,'.');
% return

%mcmc.inference='amc';
%mcmc.inference='vl';
%mcmc.inference='langevin';

%mcmc.inference='multi-amc';
%mcmc.J=512;
%mcmc.maxits=128;

mcmc.inference='ais';
mcmc.anneal='geometric';
mcmc.prop='lmc';
mcmc.nprop=1;
mcmc.J=512;
mcmc.maxits=32;

mcmc.verbose=0;

post = spm_mci_post (mcmc,M,U,Y,w_true);
disp('Posterior Mean from MCI');
disp(post.Ep)
Y_post = mci_linsqr_gen (post.Ep,M,U);

if strcmp(mcmc.inference,'ais')
    figure;plot(post.beta,mean(post.acc)); 
    title('Mean acceptance rate over trajectories');
    xlabel('Inverse temperature, \beta');
    ylabel('<Accept>');
    figure;plot(post.q);title('Normalised importance weights');
    xlabel('Sample,i');ylabel('q_i');
end

if ~strcmp(mcmc.inference,'vl')
    figure;plot(post.P(1,:),post.P(2,:),'.'); title('Posterior samples');
    xlabel('P(1)');ylabel('P(2)');
end

figure;
lw=2;
plot(M.t,Y_true,'k','LineWidth',lw);
hold on
plot(M.t,Y,'.');
plot(M.t,Y_post,'r','LineWidth',lw);

plot_post=1;
if plot_post
    % Plot posterior surface
    S.Nbins=50;
    r=2;
    S.pxy(1,:)=linspace(-r*w_true(1),r*w_true(1),S.Nbins);
    S.pxy(2,:)=linspace(-r*w_true(2),r*w_true(2),S.Nbins);
    S.param{1}='P(1)';
    S.param{2}='P(2)';
    S.name={'w_1','w_2'};
    [L,S] = mci_plot_surface (w_true,M,U,Y,S,'post');
    hold on
    ms=10;
    plot(M.pE(1),M.pE(2),'wo','MarkerSize',ms);
    plot(post.Ep(1),post.Ep(2),'wx','MarkerSize',ms);
    if ~strcmp(mcmc.inference,'vl')
        plot(post.P(1,:),post.P(2,:),'w.'); 
    end
    %plot(w_true(1),w_true(2),'k+','MarkerSize',ms);
end

spm_mci_diag(post);

stats = spm_mci_mvnpost (post,'ESS')
stats = spm_mci_mvnpost (post,'thinning')
