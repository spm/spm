function [Ep,SDp] = spm_mci_report (P,mcmc,true_P)
% Report on posterior density from MCI
% FUNCTION [Ep,SDp] = spm_mci_report (P,mcmc,true_P)
%
% P         Samples
% mcmc      Sampling options
% 
% Ep        Posterior mean
% SDp       Posterior SD
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_report.m 6548 2015-09-11 12:39:47Z will $

if nargin < 3 | isempty(true_P)
    tp=0;
else
    tp=1;
end

Np=size(P{1}.theta,1);

epoch={'Scale','Tune','Sample'};
index(1).ind=[1:mcmc.nscale];
index(2).ind=[mcmc.nscale+1:mcmc.nscale+mcmc.ntune];
index(3).ind=[mcmc.nscale+mcmc.ntune+1:mcmc.nscale+mcmc.ntune+mcmc.nsamp];

if Np==2
    figure
    for i=1:3,
        subplot(2,2,i);
        
        theta1=P{1}.theta(1,index(i).ind);
        theta2=P{1}.theta(2,index(i).ind);
        
        plot(theta1,theta2,'k.');
        grid on
        title(sprintf('%s',epoch{i}));
        xlabel('P_1');
        ylabel('P_2');
        
        hold on
        if isfield(mcmc,'init')
            Pi=mcmc.init{1};
            if i==1
                plot(Pi(1),Pi(2),'bo');
            end
        end
        if tp
            plot(true_P(1),true_P(2),'ro');
        end
    end
end

disp(sprintf('Scaling acceptance rate = %1.2f',P{1}.ar_scale));
disp(sprintf('Tuning acceptance rate = %1.2f',P{1}.ar_tune));
disp(sprintf('Sampling acceptance rate = %1.2f',P{1}.ar_sample));

% disp('Scaled covariance - used in tuning stage');
% P{1}.C
% 
% if isfield(P{1},'Ct')
%     disp('Tuned covariance - used in sampling stage');
%     P{1}.Ct
% end

% figure
% for i=1:Np,
%     subplot(Np,1,i);
%     plot(P{1}.theta(i,:));
%     grid on
%     ylabel(sprintf('P(%d)',i));
% end
% xlabel('Samples');
    
x=P{1}.theta(:,index(3).ind);
Ep=mean(x');
SDp=std(x');

%disp('Posterior mean:');
%Ep
%disp('Posterior SD:');
%SDp
