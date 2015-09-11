function [] = spm_mci_progress (x,E,i)
% Plot trajectories of parameters and neg log joint 
% FORMAT [] = spm_mci_progress (x,E,i)
% 
% x     parameters
% E     Energy = Neg Log Joint
% i     iteration number
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_progress.m 6548 2015-09-11 12:39:47Z will $

sample_win=64;
min_it=max(i-sample_win+1,1);
its=[min_it:i-1];
display(['It: ', num2str(i)]);
subplot(1,2,1); plot(its,x(its,:));
xlabel('Sample','FontName','Arial','FontSize',16);
ylabel('Params','FontName','Arial','FontSize',16);box off;
set(gca,'FontSize',16,'FontName','Arial');axis tight;

subplot(1,2,2); plot(its,E(its));
xlabel('Sample','FontName','Arial','FontSize',16);
ylabel('E','FontName','Arial','FontSize',16);box off;
set(gca,'FontSize',16,'FontName','Arial');axis tight;
drawnow;
