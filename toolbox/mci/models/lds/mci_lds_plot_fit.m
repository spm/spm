function [] = mci_lds_plot_fit (MCI,lds,n,plotfit)
% Plot fit from group LDS estimation
% FORMAT [] = mci_lds_plot_fit (MCI,lds,n,plotfit)
%
% MCI       MCI-MFX data structure
% lds       true model data structure with fields:
%
% .pinit    true init params
% .pflow    true flow params
% n         subject number
% plotfit   1 to plot model fit, 0 otherwise (default)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_lds_plot_fit.m 6697 2016-01-27 14:57:28Z spm $

if nargin < 4
    plotfit=0;
end

M=MCI.M; U=MCI.U; Y=MCI.Y;

d=size(MCI.M{1}.x0,1);

% True model responses and data
switch lds.init_par,
    case 'fixed',
        R0 = lds.pinit;
    case 'random',
        R0 = lds.pinit(:,n);
    otherwise
        % Assume known
        R0 = M{n}.x0;
end
if strcmp(lds.flow_par,'fixed')
    Pf = lds.pflow;
else
    Pf = lds.pflow(:,n);
end
M{n}.x0=R0;
y = spm_mci_fwd (Pf,M{n},U{n});

% Model fits
if plotfit
    switch MCI.assign.init_par,
        case 'fixed',
            R0 = MCI.pinit(:);
        case 'random',
            R0 = MCI.pinit_sub(:,n);
        otherwise
            % Assume known
            R0 = M{n}.x0;
    end
    if strcmp(MCI.assign.flow_par,'fixed')
        Pf = MCI.pflow(:);
    else
        Pf = MCI.pflow_sub(:,n);
    end
    M{n}.x0=R0;
    yhat = spm_mci_fwd (Pf,M{n},U{n});
end


% Data
try, ind=Y{n}.ind; catch, ind=1:M{n}.N; end
Ny=size(y,2);
ry=ceil(sqrt(Ny));

% Plotting
lw=2;
h=figure;
set(h,'Name',sprintf('Subject %d',n));
for j=1:Ny,
    subplot(ry,ry,j);
    plot(M{n}.t,y(:,j),'k','LineWidth',lw);
    hold on
    plot(M{n}.t(ind),Y{n}.y(:,j),'kx','LineWidth',lw,'MarkerSize',12);
    if plotfit
        plot(M{n}.t,yhat(:,j),'r','LineWidth',lw);
    end
    %ylim([0 5]);
    grid on
    %     if j < Ny-ry
    %         set(gca,'XTickLabel',[]);
    %     end
    
    set(gca,'FontSize',16);
    ylabel(['x_',int2str(j)]);
    if j > 2
        xlabel('Month');
    end
end
