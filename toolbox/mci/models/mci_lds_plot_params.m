function [rmse] = mci_lds_plot_params (MCI,lds)
% Plot results of group LDS estimation
% FORMAT [rmse] = mci_lds_plot_params (MCI,lds)
%
% MCI      MCI-MFX data structure
% lds      true model data structure with fields:
%
% .pinit    true init params
% .pflow    true flow params
%
% rmse      root mean square errors
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_lds_plot_params.m 6277 2014-12-04 12:16:52Z guillaume $

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
        plot(pinit(j,:),MCI.pinitK(j,:),'kx','MarkerSize',10);
        e0=MCI.pinit0-pinit;
        e2=MCI.pinitK-pinit;
    else
        hold on
        plot(pinit(j,:),MCI.pinitK(j),'kx','MarkerSize',10);
        e0=MCI.pinit0(:)*ones(1,Nsub)-pinit;
        e2=MCI.pinitK(:)*ones(1,Nsub)-pinit;
    end
    plot([min_w max_w],[min_w max_w],'k-','LineWidth',2);
    set(gca,'FontSize',18);
    xlabel('True');
    ylabel('Estimated');
    grid on
    
end

if use_median
    rmse.pinit0=sqrt(median(diag(e1'*e1))/d);
    rmse.pinitK=sqrt(median(diag(e2'*e2))/d);
else
    rmse.pinit0=sqrt(mean(diag(e0'*e0))/d);
    rmse.pinitK=sqrt(mean(diag(e2'*e2))/d);
end

fprintf('\nInitial state parameters:\n');
fprintf('Pinit0 RMSE=%1.2f\n',rmse.pinit0);
fprintf('PinitK RMSE=%1.2f\n',rmse.pinitK);
  


% Flow estimates
h=figure;
set(h,'Name','Flow Parameters');
min_v=min(pflow);
max_v=max(pflow);
plot([min_v max_v],[min_v max_v],'k-','LineWidth',2);
hold on
if strcmp(MCI.assign.flow_par,'random')
    if strcmp(lds.flow_par,'fixed')
        e2=MCI.pflowK-pflow*ones(1,Nsub);
        for n=1:Nsub,
            plot(pflow,MCI.pflowK(:,n),'kx','MarkerSize',10);
        end
    else
        e2=MCI.pflowK-pflow;
        for n=1:Nsub,
            plot(pflow(:,n),MCI.pflowK(:,n),'kx','MarkerSize',10);
        end
    end
else
    plot(pflow,MCI.pflowK,'kx','MarkerSize',10);
    if strcmp(lds.flow_par,'fixed')
        e2=MCI.pflowK-pflow';
    else
        e2=MCI.pflowK-pflow*ones(1,Nsub);
    end
end
set(gca,'FontSize',18);
xlabel('True');
ylabel('Estimated');

Np=length(pflow);
if use_median
    rmse.pflowK=sqrt(median(diag(e2'*e2))/Np);
else
    rmse.pflowK=sqrt(mean(diag(e2'*e2))/Np);
end

fprintf('\nFlow parameters:\n');
fprintf('PflowK RMSE=%1.2f\n',rmse.pflowK);
