
clear all
close all

disp('Synthetic data from nonlinear ODE model');
disp('defined in Ramsay et al. (2007)');
disp('based on Van der Pol oscillator and which');
disp('reduces to Fitzhugh-Nagumo for certain parameters');

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

tic;
Y = mci_ramsay_gen(P,M,U);
toc

mci_plot_outputs(M,Y);

%mcmc.inference='vl';
%mcmc.inference='langevin';
%mcmc.maxits=16;
%mcmc.verbose=1;

mcmc.inference='ais';
mcmc.anneal='geometric';
mcmc.prop='lmc';
mcmc.nprop=1;
mcmc.J=16;
mcmc.maxits=8;
mcmc.rec_traj=1;

post = spm_mci_post (mcmc,M,U,Y,P);

if ~strcmp(mcmc.inference,'vl')
    spm_mci_diag(post);
    
    stats = spm_mci_mvnpost (post,'ESS')
    stats = spm_mci_mvnpost (post,'thinning')
    
    for j=1:length(P),
        spm_mci_quantiles (post,j,0);
    end
end
    
load ramsay-surface
figure
surf(S.x,S.y,L);
xlabel(S.name{1});
ylabel(S.name{2});

figure;
imagesc(S.pxy(1,:),S.pxy(2,:),L);
axis xy
hold on
xlabel(S.name{1});
ylabel(S.name{2});
hold on
ms=10;
plot(M.pE(1),M.pE(2),'wo','MarkerSize',ms);
plot(P(1),P(2),'w+','MarkerSize',ms);
if ~strcmp(mcmc.inference,'vl')
    j=post.ind;
    plot(post.P(1,j),post.P(2,j),'w.');
end
plot(post.Ep(1),post.Ep(2),'kx','MarkerSize',ms);

if strcmp(mcmc.inference,'ais')
    % plot individual trajectories
    figure
    xlim([-2 0]);
    ylim([-2 0]);
    hold on
    for i=1:mcmc.maxits,
        xx=squeeze(post.traj(i,1,:));
        yy=squeeze(post.traj(i,2,:));
        plot(xx,yy,'k-');
        plot(xx(end),yy(end),'kx','MarkerSize',ms);
        xlabel(S.name{1});
        ylabel(S.name{2});
        disp('Press space to see more trajectories');
        pause
        
    end
end
