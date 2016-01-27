function [] = mci_plot_dist_multi (dist,name,P) 
% Plot (multiple) densities
% FORMAT [] = mci_plot_dist_multi (dist,name,P) 
% 
% dist{i}   ith distribution 
%
% .Ep       mean
% .P        [Np x Ns] sample matrix
% .ind      indices of samples (eg. post burn-in)
% .names    names of variables
% .color    eg 'r','k','b'
% .order    eg. [1,3,4,2] to plot only variables 1,3,4 and 2 in that order
%
% name      name of parameters
% P         true parameters (optional)
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_plot_dist_multi.m 6697 2016-01-27 14:57:28Z spm $

Ndist=length(dist);
for i=1:Ndist
    if strcmpi(dist{i}.type,'sample')
        figure;
        plot(dist{i}.P');
        xlabel('Sample');
        ylabel('Parameters');
        title(sprintf('%s trajectories',name));
    end
end

Np=length(dist{1}.Ep);
try
    order=dist{1}.order;
    Np=length(order);
catch
    order=[1:Np];
end

% Plot univariate densities
rNp=ceil(sqrt(Np));
h=figure;
set(h,'Name',name);
for j=1:Np,
    subplot(rNp,rNp,j);
    jp=order(j);
    for i=1:Ndist,
        mci_plot_dist(dist{i},jp);
        hold on
    end
    yl=get(gca,'YLim');
    if nargin > 2
        plot([P(jp),P(jp)],[0,yl(2)],'r','LineWidth',2);
    end
end
if nargin > 2
    disp('True parameters shown in red');
end

plot_bivariate=0;
if plot_bivariate
    % Plot bivariate densities
    rNp=ceil(sqrt(Np));
    h=figure;
    set(h,'Name',name);
    k=1;
    q=dist{1}.P(:,dist{1}.ind);
    for i=1:Np,
        for j=1:Np,
            if j > i
                subplot(Np,Np,k);
                plot(q(i,:),q(j,:),'k.');
                xlabel(dist{1}.names{i});
                ylabel(dist{1}.names{j});
                hold on
                plot(P(i),P(j),'ro');
            end
            k=k+1;
        end
    end
    if nargin > 2
        disp('True parameters shown in red');
    end
    
    % Posterior Correlation Matrix
    for i=1:Ndist,
        if strcmp(lower(dist{i}.type),'sample')
            disp(' ');
            disp(sprintf('Distribution %d',i));
            disp('Posterior Correlation Matrix:');
            corrcoef(dist{i}.P(:,dist{i}.ind)')
        end
    end
end
