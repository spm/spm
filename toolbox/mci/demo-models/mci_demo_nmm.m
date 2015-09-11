
clear all
close all

% Number of parameters to estimate
Np=2;
%Np=6;
%Np=21;

disp('Synthetic data from two region Neural Mass Model');
disp(sprintf('Model has %d parameters',Np));

% Backward connection ?
back=1;

% Observation noise SD
sd=0.01;

[M,U] = mci_nmm_struct (back,sd,Np);

%M.Ce=diag([0.03^2,0.01^2]);

% Parameters
switch Np
    case 2
        disp('Estimating extrinsic connections');
        P=[1,1]';
        %P=[0.5,1]';
    case 6
        disp('Estimating intrinsic and extrinsic connections');
        P=[1,1,0,0,1,0]';
    case 21,
        disp('Estimating all connections');
        % All params set to prior mean except f/b
        P=M.pE;
        P.A{1}(2,1)=1; % Forward connection
        P.A{2}(1,2)=1; % Backward connection
end

Y = mci_nmm_gen (M,U,P);

mci_plot_outputs(M,Y);

%mcmc.inference='amc';
%mcmc.inference='vl';
mcmc.inference='langevin';

mcmc.maxits=1024;
mcmc.verbose=0;

post = spm_mci_post (mcmc,M,U,Y,P);

switch mcmc.inference
    case 'langevin',
        ind=[1:size(post.Ce,3)];
        mci_plot_noiseSD (post.Ce,ind);
    case 'vl',
        disp('Error SD:');
        disp(sqrt(diag(post.Ce)));
end

diag.traceplot=1;
diag.eplot=1;
if strcmp(mcmc.inference,'langevin')
    diag.bplot=1;
else
    diag.bplot=0;
end
spm_mci_diag(post,diag);

disp(sprintf('Elapsed time:%1.2f',post.els));

% Plot posterior surface
plot_post=0;
if plot_post
    S.Nbins=20;
    S.pxy(1,:)=linspace(-0.2,1.5,S.Nbins);
    S.pxy(2,:)=linspace(-0.2,1.5,S.Nbins);
    S.param{1}='P(1)';
    S.param{2}='P(2)';
    S.name={'P(1)','P(2)'};
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
end

stats = spm_mci_mvnpost (post,'ESS')
stats = spm_mci_mvnpost (post,'thinning')

spm_mci_quantiles (post,1,0);
xlabel('\theta_1');
ylabel('p(\theta_1|Y)');

spm_mci_quantiles (post,2,0);
xlabel('\theta_2');
ylabel('p(\theta_2|Y)');

%[x,pnum,pgauss]=spm_mci_postslices (post,M,U,Y);
