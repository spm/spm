function [y] = spm_mci_quantiles (post,j,q3,expP)
% Plot histogram and quantiles of posterior density
% FORMAT [y] = spm_mci_quantiles (post,j,q3,expP)
%
% post      posterior data structure
% j         jth variate
% q3        plot quantiles on histogram
% expP      exponentiate parameters before plotting ?
%
% y         2.5%, 50%, 97.5% quantiles
%
% Solid lines show quantiles from posterior samples
% Dotted lines under Gaussian assumptions
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_quantiles.m 6697 2016-01-27 14:57:28Z spm $

try, plotq=q3; catch, plotq=1; end
try, expP=expP; catch, expP=0; end

q = [.025 .5 .975];
N = length(q);
x = post.P(j,post.ind);
if expP, x=exp(x); end
y = quantile(x,q);
%figure
hist(x);
hold on
yl=get(gca,'YLim');
ym=yl(2);

h = findobj(gca,'Type','patch');
set(h,'FaceColor', [.5 .5 .5]);

m=mean(x);
s=std(x);
yg=spm_invNcdf(q,m,s^2);

%cols={'k','b','r','b','k'};
if plotq
    lw=2;
    for i = 1:N,
        plot([y(i) y(i)],[0 ym],'k','LineWidth',lw);
        plot([yg(i) yg(i)],[0 ym],'k:','LineWidth',lw);
    end
end
