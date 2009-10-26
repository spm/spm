function out=spm_api_bmc(F,N,exp_r,xp,family)
% API to select and compare DCMs using Bayesian model comparison
% FORMAT out=spm_api_bmc(F,N,alpha,exp_r,xp)
%
% INPUT:
% F      - Matrix/Vector of log model evidences
% N      - vector of model names
% alpha  - vector of model probabilities
% exp_r  - expectation of the posterior p(r|y)
% xp     - exceedance probabilities
%
% OUTPUT:
% out    - conditional probability of DCMs (when using fixed effect method)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_api_bmc.m 3508 2009-10-26 13:04:17Z maria $

if isempty(xp)
    inf_method = 'FFX';
else
    inf_method = 'RFX';
end

if nargin>4
    plot_family = 1;
else
    plot_family = 0;
end

nm = length(N);

Fgraph  = spm_figure('GetWin','Graphics','BMS: results');

switch inf_method

    %======================================================================
    % Fixed Effect
    %======================================================================
    case ('FFX')
        
        figure(Fgraph);
        %-Compute conditional probability of DCMs under flat priors.
        %------------------------------------------------------------------
        F    = F - min(F);
        i    = F < (max(F) - 32);
        P    = F;
        P(i) = max(F) - 32;
        P    = P - min(P);
        P    = exp(P);
        P    = P/sum(P);

        %-Display results
        %------------------------------------------------------------------
        subplot(2,1,1)
        bar(1:nm,F)
        set(gca,'XTick',1:nm)
        set(gca,'XTickLabel',1:nm)
        ylabel('Log-evidence (relative)','Fontsize',14)
        xlabel('Models','Fontsize',14)
        title({'Bayesian Model Selection'},'Fontsize',14)
        axis square
        grid on
        
        subplot(2,1,2)
        bar(1:nm,P)
        set(gca,'XTick',1:nm)
        set(gca,'XTickLabel',1:nm)
        ylabel('Posterior Model Probability','Fontsize',14)
        title('Bayesian Model Selection','Fontsize',14)
        xlabel('Models','Fontsize',14)
        axis square
        grid on
        
        out = P;
        
        if plot_family
            
            Ffam = spm_figure('CreateWin','Graphics','BMS: results');
            figure(Ffam);
           
            %-Display results - families
            %------------------------------------------------------------------
            subplot(2,1,1)
            Nfam = length(family.post);
            bar(1:nm,P)
            set(gca,'XTick',1:nm)
            set(gca,'XTickLabel',family.names)
            ylabel('Posterior Model Probability','Fontsize',14)
            title('Bayesian Model Selection','Fontsize',14)
            xlabel('Models','Fontsize',14)
            axis square
            grid on
            
            subplot(2,1,2)
            bar(1:Nfam,family.post)
            set(gca,'XTick',1:Nfam)
            set(gca,'XTickLabel',family.names)
            ylabel('Posterior Family Probability','Fontsize',14)
            xlabel('Families','Fontsize',14)
            title('Bayesian Model Selection','Fontsize',14)
            axis square
            grid on
            
        end
    
    %======================================================================
    % Random Effect
    %======================================================================
    case ('RFX')

        %-Display results
        %------------------------------------------------------------------
        subplot(2,1,1)
        bar(1:length(N),exp_r)
        set(gca,'XTick',1:length(N))
        set(gca,'XTickLabel',1:nm)
        ylabel('Expected Posterior Probability','Fontsize',14)
        xlabel('Models','Fontsize',14)
        title('Bayesian Model Selection','Fontsize',14)
        axis square
        grid on

        subplot(2,1,2)
        bar(1:length(N),xp')
        set(gca,'XTick',1:length(N))
        set(gca,'XTickLabel',1:nm)
        ylabel('Exceedance Probability','Fontsize',14)
        xlabel('Models','Fontsize',14)
        title('Bayesian Model Selection','Fontsize',14)
        axis square
        grid on
        
        out = [];
        
        if plot_family
            
            Ffam = spm_figure('CreateWin','Graphics','BMS: results');
            figure(Ffam);
            
            %-Display results - families
            %------------------------------------------------------------------
            subplot(2,1,1)
            Nfam = length(family.exp_r);
            bar(1:Nfam,family.exp_r)
            set(gca,'XTick',1:Nfam)
            set(gca,'XTickLabel',family.names)
            ylabel('Posterior Family Probability','Fontsize',14)
            xlabel('Families','Fontsize',14)
            title('Bayesian Model Selection','Fontsize',14)
            axis square
            grid on
            
            subplot(2,1,2)
            bar(1:Nfam,family.xp')
            set(gca,'XTick',1:Nfam)
            set(gca,'XTickLabel',family.names)
            ylabel('Exceedance Family Probability','Fontsize',14)
            xlabel('Families','Fontsize',14)
            title('Bayesian Model Selection','Fontsize',14)
            axis square
            grid on
            
        end
        
end


