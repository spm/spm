function out=spm_api_bmc(F,N,alpha,exp_r,xp)
% API to select and compare DCMs using Bayesian model comparison
% FORMAT out=spm_api_bmc(F,N,alpha,exp_r,xp)
%
%INPUT :
% F     - Matrix/Vector of log model evidnces
%
%
%__________________________________________________________________________
% 
%Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_api_bmc.m 2628 2009-01-20 16:32:17Z maria $

if nargin < 3
    inf_method='FFX';
else
    inf_method='RFX';
end
nm = length(N);


switch inf_method

    case ('FFX')
        % compute conditional probability of DCMs under flat priors.
        %--------------------------------------------------------------------------
        F    = F - min(F);
        i    = F < (max(F) - 32);
        P    = F;
        P(i) = max(F) - 32;
        P    = P - min(P);
        P    = exp(P);
        P    = P/sum(P);

        % display results
        %--------------------------------------------------------------------------
        Fgraph  = spm_figure('GetWin','Graphics'); figure(Fgraph); clf

        
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
        out=P;
        
    case ('RFX')
        
          Fgraph  = spm_figure('GetWin','Graphics'); figure(Fgraph); clf

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
            out=[];
        
end
