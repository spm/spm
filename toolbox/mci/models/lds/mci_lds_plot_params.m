function [rmse] = mci_lds_plot_results (MCI,lds)
% Plot results of group LDS estimation
% FORMAT [rmse] = mci_lds_plot_results (MCI,lds)
%
% MCI      MCI-MFX data structure
% lds      true model data structure with fields:
%
% .pinit    true init params
% .pflow    true flow params
%
% rmse      root mean square errors
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_lds_plot_params.m 6548 2015-09-11 12:39:47Z will $

pinit=lds.pinit; pflow=lds.pflow;
d=size(MCI.M{1}.x0,1);
Nsub=length(MCI.Y);

% Use mean or median when computing RMSE
use_median=0;

% Initial state estimates
for j=1:d,
    h=figure;
    set(h,'Name',sprintf('Initial State %d',j));
    
    min_w=min(pinit(j,:));
    max_w=max(pinit(j,:));
    
    if strcmp(MCI.assign.init_par,'random')
        hold on
        plot(pinit(j,:),MCI.pinit_sub(j,:),'kx','MarkerSize',10);
        e2=MCI.pinit_sub-pinit;
    else
        hold on
        plot(pinit(j,:),MCI.pinit(j),'kx','MarkerSize',10);
        e2=MCI.pinit(:)*ones(1,Nsub)-pinit;
    end
    plot([min_w max_w],[min_w max_w],'k-','LineWidth',2);
    set(gca,'FontSize',18);
    xlabel('True');
    ylabel('Estimated');
    grid on
    
end

if use_median
    rmse.pinit=sqrt(median(diag(e2'*e2))/d);
else
    rmse.pinit=sqrt(mean(diag(e2'*e2))/d);
end

disp(' ');
disp('Initial state parameters:');
disp(sprintf('Pinit RMSE=%1.2f',rmse.pinit));
  


% Flow estimates
h=figure;
set(h,'Name','Flow Parameters');
min_v=min(pflow);
max_v=max(pflow);
plot([min_v max_v],[min_v max_v],'k-','LineWidth',2);
hold on
if strcmp(MCI.assign.flow_par,'random')
    if strcmp(lds.flow_par,'fixed')
        e2=MCI.pflow-pflow*ones(1,Nsub);
        for n=1:Nsub,
            plot(pflow,MCI.pflow(:,n),'kx','MarkerSize',10);
        end
    else
        e2=MCI.pflowK-pflow;
        for n=1:Nsub,
            plot(pflow(:,n),MCI.pflow(:,n),'kx','MarkerSize',10);
        end
    end
else
    plot(pflow,MCI.pflow,'kx','MarkerSize',10);
    if strcmp(lds.flow_par,'fixed')
        e2=MCI.pflow-pflow;
    else
        e2=MCI.pflow-pflow*ones(1,Nsub);
    end
end
set(gca,'FontSize',18);
xlabel('True');
ylabel('Estimated');

Np=length(pflow);
if use_median
    rmse.pflow=sqrt(median(diag(e2'*e2))/Np);
else
    rmse.pflow=sqrt(mean(diag(e2'*e2))/Np);
end

disp(' ');
disp('Flow parameters:');
disp(sprintf('Pflow RMSE=%1.2f',rmse.pflow));