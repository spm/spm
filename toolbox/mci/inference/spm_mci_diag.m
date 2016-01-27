function [mess] = spm_mci_diag (post,diag)
% Monte Carlo Diagnostics
% FORMAT [mess] = spm_mci_diag (post,diag)
%
% post      posterior distribution
% diag      diagnostic info
%           .ind            indices of samples to look at
%           .traceplot      (1/0) for trace plots
%           .autoplot       (1/0) for autocorrelations
%           .essplot        (1/0) for effective sample sizes
%           .eplot          (1/0) for energy (neg log joint) traj
%           .bplot          (1/0) for Bayes factor of f/b transitions
%
% ess      effective sample size (for each parameter)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_diag.m 6697 2016-01-27 14:57:28Z spm $

Ns=size(post.P,2);
try, traceplot=diag.traceplot; catch, traceplot=1; end
try, autoplot=diag.autoplot; catch, autoplot=0; end
try, essplot=diag.essplot; catch, essplot=0; end
try, eplot=diag.eplot; catch, eplot=1; end
try, bplot=diag.bplot; catch, bplot=1; end
try, ind=diag.ind; catch, ind=[1:Ns]; end

P=post.P(:,ind);

[Np,Ns]=size(P);
rp=ceil(sqrt(Np));

if eplot
    h=figure;
    set(h,'Name','Energy Trajectory');
    %plot(post.Etraj(ind));
    set(gca,'FontSize',16);
    semilogx(post.E(ind),'k');
    grid on
    xlabel('Sample');
    ylabel('Energy');
end

if bplot
    j=ind(2:end);
    
    h=figure;
    set(h,'Name','Proposal Details');
    %plot(post.bayes_fb(ind));
    
    if isfield(post,'bayes_fb')
        subplot(2,2,1);
        semilogx(post.bayes_fb(j),'k');
        set(gca,'FontSize',16);
        grid on
        xlabel('Sample');
        ylabel('LogBF (FB)');
    end
    
    if isfield(post,'dL')
        subplot(2,2,3);
        semilogx(post.dL(j),'k');
        set(gca,'FontSize',16);
        grid on
        xlabel('Sample');
        ylabel('dL');
    end
    
    if isfield(post,'dL') && isfield(post,'bayes_fb')
        subplot(2,2,2);
        semilogx(post.dL(j)-post.bayes_fb(j),'k');
        set(gca,'FontSize',16);
        grid on
        xlabel('Sample');
        ylabel('Log r');
    end
    
    if isfield(post,'acc')
        subplot(2,2,4);
        semilogx(post.acc(j),'k');
        ylim([-0.1 1.1]);
        set(gca,'FontSize',16);
        grid on
        xlabel('Sample');
        ylabel('Accepted');
    end
end

if traceplot
    % Trace plots for each parameter
    h=figure;
    set(h,'Name','Parameter Samples');
    for p=1:Np,
        subplot(rp,rp,p);
        plot(P(p,:),'k');
        %set(gca,'FontSize',16);
        xlabel('Sample');
        ylabel(sprintf('P(%d)',p));
        grid on
    end
end

if autoplot
    % Autocorrelation plots for each parameter
    h=figure;
    set(gca,'FontSize',16);
    set(h,'Name','Parameter Autocorrelations');
    Nlags=ceil(Ns/10);
    for p=1:Np,
        subplot(rp,rp,p);
        plot(autocorr(P(p,:),Nlags),'k');
        xlabel('Sample');
        ylabel(sprintf('P(%d)',p));
        grid on
    end
end

if essplot
    % Effective sample sizes
    h=figure;
    set(h,'Name','Effective Sample Size');
    for p=1:Np,
        ess(p)=spm_mci_ess(P(p,:));
    end
    mess=mean(ess);
    disp('Mean Effective Sample Size:');
    disp(mess);
    
    bar(ess);
    xlabel('Parameter');
    ylabel('ESS');
else
    mess=[];
end

