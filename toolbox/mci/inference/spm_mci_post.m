function [post] = spm_mci_post (inference,M,U,Y,true_P,verbose)
% Estimate posterior density
% FORMAT [post] = spm_mci_post (inference,M,U,Y,true_P,verbose)
%
% inference    'mh','vl','langevin' 
% M             model structure
% U             inputs (shouldn't be empty)
% Y             data
% true_P        true parameters (if known)
% verbose       0/1 to plot progress (default 0)
%
% post          structure containing posterior (mean, samples etc)
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_post.m 6275 2014-12-01 08:41:18Z will $

try tp = true_P; catch tp=0; end
try verbose=verbose; catch verbose=0; end

tic;
switch inference,
    case 'mh',
        disp('Estimating parameters using Metropolis-Hastings');
        disp(' ');
        
        % MH defaults
        Nsamp=200;
        mcmc = spm_mci_popdef (Nsamp,Nsamp,Nsamp);
        mcmc.verbose=verbose;
        
        % Draw initialisation point from prior ?
        %mcmc.init{1}=spm_normrnd(M.pE,M.pC,1);
        mcmc.init{1}=M.pE;
        
        MM{1}=M;
        UU{1}=U;
        tic;
        [Psamp,logev,D,MM] = spm_mci_pop (mcmc,MM,UU,Y);
        toc
        M=MM{1};
        
        if tp
            [post.Ep,post.SDp]=spm_mci_report (Psamp,mcmc,true_P);
            disp('True params:');
            disp(true_P)
        else
            disp('Initial params:');
            disp(mcmc.init{1});
            [post.Ep,post.SDp]=spm_mci_report (Psamp,mcmc);
        end
        
        post.Ep=post.Ep';
        post.Psamp=Psamp;
        post.logev=logev;
        post.D=D;
        post.mcmc=mcmc;
        
    case 'vl',
        disp('Estimating parameters using Variational Laplace');
        disp(' ');
        
        D.y=Y;
        if ~verbose
            M.nograph=1;
        end
        
        % Check to see if VL will run
        UI.u=U';
        UI.dt=M.T/M.N;
        %tmp=spm_int(M.pE,M,UI);

        M.hE=9;
        M.hC=0.1;
        %M.P = spm_normrnd(M.pE,M.pC,1);
        [Ep,Cp,Ch] = spm_nlsi_GN (M,UI,D);
        
        %post.Ep=M.V'*spm_vec(Ep);
        post.Ep=spm_vec(Ep);
        post.Cp=Cp;
        
        % Eigenrep needed for spm_mci_joint later
        M = spm_mci_minit (M);
        
    case 'langevin',
        disp('Estimating parameters using Langevin Monte Carlo');
        disp(' ');
        
        mcmc.maxits=512;
        mcmc.verbose=verbose;
        [M,stats]=spm_mci_lgv(mcmc,M,U,Y);
        Psamp=stats.P';
        Nsamp=size(Psamp,1);
        burn_in=round(0.3*size(Psamp,1));
        post.ind=[burn_in+1:Nsamp];
        post.targ=Psamp(post.ind,:);
        post.Ep=mean(post.targ)';
        post.Cp=cov(post.targ);
        post.P=stats.P;
        post.Ce=stats.Ce;
        post.Etraj=stats.E;
        post.dEdit=stats.dEdit;
        
    otherwise
        disp('Unknown inference method');
end

post.els=toc;
disp(sprintf('Optimisation time = %1.2f seconds',post.els));

% Generate predictions
if isfield(M,'IS')
    post.Yhat = feval(M.IS,post.Ep,M,U);
else
    post.Yhat = spm_mci_fwd (post.Ep,M,U);
end

figure
rm=ceil(sqrt(M.l));
lw=2;
for i=1:M.l,
    if M.l>3
        subplot(rm,rm,i);
    else
        subplot(M.l,1,i);
    end
    plot(M.t,Y(:,i),'LineWidth',lw);
    hold on
    plot(M.t,post.Yhat(:,i),'r','LineWidth',lw);
    grid on
    set(gca,'FontSize',16);
    legend('Data','Fit');
    xlabel('Time');
    ylabel(sprintf('y(%d)',i));
end

% get parameters in reduced space
Pr=M.V'*(post.Ep-M.vpE);
post.L_est = spm_mci_joint (Pr,M,U,Y);
disp(sprintf('Estimated Log Joint=%1.2f',post.L_est));

Pr=M.V'*(post.Ep-M.vpE);
pt=4;
if M.Np > pt
    hp=figure;
    set(hp,'Name','Parameters');
    plot(Pr,'r','LineWidth',lw);
    xlabel('Parameter');
    set(gca,'FontSize',16);
    grid on
else
    disp('Estimated (latent) params:');
    disp(Pr);
end

if tp
    % get parameters in reduced space
    Pr=M.V'*(spm_vec(true_P)-M.vpE);
    
    if M.Np > pt
        hold on
        plot(Pr,'LineWidth',lw);
        legend('Estimated','True');
    else
        disp('True (latent) params:');
        disp(Pr);
    end
    
    post.L_true = spm_mci_joint (Pr,M,U,Y);
    disp(sprintf('True Log Joint=%1.2f',post.L_true));
    disp(' ');
end

switch inference,
    case {'mh','langevin'},
        post.type='sample';
    otherwise
        post.type='gaussian';
end