function [] = mci_plot_dist (dist,j,xlims) 
% Plot probability density
% FORMAT [] = mci_plot_dist (dist,j,xlims) 
% 
% dist      struct with fields
%
% .Ep       posterior mean
% .P        [Np x Ns] sample matrix
% .ind      indices of samples dist burn-in
% .names
% .ks       set to 1 for kernel smoothing (default)
% j         jth variable
% xlims     xlims(1,2) for lower/upper limits
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: mci_plot_dist.m 6697 2016-01-27 14:57:28Z spm $

try, ks=dist.ks; catch, ks=1; end

if nargin > 2
    limits=1;
else
    limits=0;
end
lw=2;

switch lower(dist.type),
    case 'boxplot',
        boxplot(dist.P(j,dist.ind));
    case 'sample',
        if ks
            [g,xi]=ksdensity(dist.P(j,dist.ind));
            plot(xi,g,dist.color,'LineWidth',lw);set(gca,'YTick',[]);
        else
            hist(dist.P(j,dist.ind),20);
        end
    case 'gaussian',
        m=dist.Ep(j);s=sqrt(full(dist.Cp(j,j)));
        xi=linspace(m-4*s,m+4*s,100);
        g=spm_Npdf(xi,m,s^2);
        plot(xi,g,dist.color,'LineWidth',lw);set(gca,'YTick',[]);
end
set(gca,'FontSize',16);
xlabel(dist.names{j});
grid on

if limits
    xlim([xlims(1) xlims(2)]);
end