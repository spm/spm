function [post] = spm_mci_post (mcmc,M,U,Y,true_P)
% Estimate posterior density
% FORMAT [post] = spm_mci_post (mcmc,M,U,Y,true_P)
%
% mcmc          .inference = 'amc','ais','vl' or 'langevin' 
%               .verbose = 0 or 1 to plot progress (default 0)
%               .maxits = max number of iterations for sampling
%               .init = init parameter values (default is prior mean)
% M             model structure
% U             inputs (shouldn't be empty)
% Y             data
% true_P        true parameters (if known)
%
% post          structure containing posterior (mean, samples etc)
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_post.m 6697 2016-01-27 14:57:28Z spm $

try, verbose=mcmc.verbose; catch, verbose=0; end
if nargin < 5 || isempty(true_P)
    tp=0;
else
    tp=1;
end

tic;
switch mcmc.inference,
    
    case 'ais',
        disp('Annealed Importance Sampling');
        disp(' ');
        
        mcmc.nsamp=mcmc.maxits;
        
        tic;
        post = spm_mci_ais (mcmc,M,U,Y);
        toc
        
        Nsamp=size(post.P,2);
        post.ind=[1:Nsamp];
        post.Ep=mean(post.P(:,post.ind)')';
        post.mcmc=mcmc;

        % Eigenrep needed for spm_mci_joint later
        M = spm_mci_minit (M);
        
    case 'multi-amc',
        disp('Multiple chains of Adaptive Monte Carlo');
        disp(' ');
        Nsamp=ceil(0.5*mcmc.J);
        Nscale=0;
        Ntune=Nsamp;
        mc = spm_mci_popdef (Nscale,Ntune,Nsamp);
        mc.verbose=verbose;
        
        MM{1}=M;
        UU{1}=U;
        
        for it=1:mcmc.maxits,
            % Loop over chains
            mcmc.init{1}=spm_normrnd(M.pE,M.pC,1);
            Psamp = spm_mci_pop (mc,MM,UU,Y);
            if it ==1
                P=Psamp{1}.theta;
            else
                P=[P,Psamp{1}.theta];
            end
        end
        post.ind=[1:size(P,2)];
        post.Ep=mean(P(:,post.ind)')';
        post.P=P;
        post.mcmc=mcmc;
        
        % Eigenrep needed for spm_mci_joint later
        M = spm_mci_minit (M);
        
    case 'amc',
        disp('Adaptive Monte Carlo');
        disp(' ');
            
        % MH defaults
        tune=1;
        if tune
            Nscale=ceil(0.25*mcmc.maxits);
            Ntune=Nscale;
            Nsamp=ceil(0.5*mcmc.maxits);
        else
            Nscale=ceil(0.5*mcmc.maxits);
            Ntune=0;
            Nsamp=Nscale;
        end
        mc = spm_mci_popdef (Nscale,Ntune,Nsamp);
        mc.verbose=verbose;
        
        % Draw initialisation point from prior ?
        %mcmc.init{1}=spm_normrnd(M.pE,M.pC,1);
        try, mc.init{1}=mcmc.init; catch, mc.init{1}=spm_vec(M.pE); end
        
        MM{1}=M;
        UU{1}=U;
        tic;
        [Psamp,logev,D,MM] = spm_mci_pop (mc,MM,UU,Y);
        toc
        M=MM{1};
        
        if tp
            [post.Ep,post.SDp]=spm_mci_report (Psamp,mc,true_P);
        else
            disp('Initial params:');
            disp(mc.init{1});
            [post.Ep,post.SDp]=spm_mci_report (Psamp,mc);
        end
        
        post.Ep=post.Ep';
        post.P=Psamp{1}.theta;
        post.E=-Psamp{1}.logq;
        post.dL=Psamp{1}.dL;
        post.bayes_fb=zeros(1,length(post.E));
        post.acc=Psamp{1}.acc;
        post.logev=logev;
        post.D=D;
        post.mcmc=mcmc;
        post.ind=mc.ind_samp;
        
    case 'vl',
        disp('Variational Laplace');
        disp(' ');
        
        D.y=Y;
        if ~verbose
            M.nograph=1;
        end
        
        if isstruct(U)
            UI=U;
        else
            UI.u=U';
            UI.dt=M.T/M.N;
        end
        
        % Change VL defaults
        if isfield(mcmc,'maxits'), M.Nmax=mcmc.maxits; end
        if isfield(mcmc,'init'), M.P=mcmc.init; end
        if isfield(mcmc,'hE'), M.hE=mcmc.hE; end
        if isfield(mcmc,'hC'), M.hC=mcmc.hC; end
        
        
        [Ep,Cp,Eh,F,tmp1,tmp2,tmp3,k] = spm_nlsi_GN (M,UI,D);
        
        post.Ep=spm_vec(Ep);
        post.Cp=Cp;
        post.Eh=Eh;
        post.Ce=diag(1./exp(Eh));
        post.logev=F;
        post.its=k;
        if isfield(M,'P')
            post.init=M.P;
        end
        
        % Eigenrep needed for spm_mci_joint later
        M = spm_mci_minit (M);

    case 'langevin',
        disp('Langevin Monte Carlo');
        disp(' ');
        
        [M,stats]=spm_mci_lgv(mcmc,M,U,Y);
        Psamp=stats.P';
        Nsamp=size(Psamp,1);
        burn_in=round(0.3*size(Psamp,1));
        post=stats;
        post.ind=[burn_in+1:Nsamp];
        post.targ=Psamp(post.ind,:);
        post.Ep=mean(post.targ)';
        
        % 2.5%, 50% and 97.5% quantiles
        q = [.025 .5 .975];
        for j=1:M.Np,
            sx = post.P(j,post.ind);
            post.quantiles(j,:) = quantile(sx,q);
        end
        
    otherwise
        disp('Unknown inference method');
end

if strcmp(mcmc.inference,'VL')
    Up=UI;
else
    Up=U;
end

% Generate data fit from posterior mean
if isfield(M,'IS')
    if strcmp(M.IS,'spm_gen_erp')
        Ps=spm_unvec(post.Ep,M.pE);
        post.Yhat=feval('spm_gen_erp',Ps,M,Up);
    else
        post.Yhat = feval(M.IS,post.Ep,M,Up);
    end
else
    post.Yhat = spm_mci_fwd (post.Ep,M,Up);
end

post.els=toc;
disp(sprintf('Optimisation time = %1.2f seconds',post.els));

lw=2;
if ~isfield(M,'IS')
    % Plot time series
    figure
    rm=ceil(sqrt(M.l));
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
end

% get parameters in reduced space
Pr=M.V'*(post.Ep-M.vpE);
post.L_est = spm_mci_joint (Pr,M,U,Y);
disp(sprintf('Estimated Log Joint=%1.2f',post.L_est));

if tp    
    % get parameters in reduced space
    Pr=M.V'*(spm_vec(true_P)-M.vpE);
    post.L_true = spm_mci_joint (Pr,M,U,Y);
    disp(sprintf('True Log Joint=%1.2f',post.L_true));
    disp(' ');
end

pt=4;
if M.Np > pt
    hp=figure;
    set(hp,'Name','Parameters');
    plot(post.Ep,'r','LineWidth',lw);
    xlabel('Parameter');
    set(gca,'FontSize',16);
    grid on
else
    disp('Estimated (latent) params:');
    disp(post.Ep);
end

if tp
    if M.Np > pt
        hold on
        plot(spm_vec(true_P),'LineWidth',lw);
        legend('Estimated','True');
    else
        disp('True (latent) params:');
        disp(spm_vec(true_P));
    end
end

switch mcmc.inference,
    case {'amc','langevin'},
        post.type='sample';
    otherwise
        post.type='gaussian';
end