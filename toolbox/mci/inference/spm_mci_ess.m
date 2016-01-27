function [ess,m] = spm_mci_ess (x,p)
% Compute Effective Sample Size
% FORMAT [ess,m] = spm_mci_ess (x,p)
%
% x      Univariate time series
% p      Maximum lag for autocovariance estimation
%
% ess    Effective Sample Size
% m      Number of lags used in ESS estimate
%
% This routine is based on the Initial Positive Sequence estimate
% proposed in C. Geyer (1992) Practical Markov Chain Monte Carlo, 
% Statistical Science, 7(4):473-511.
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mci_ess.m 6697 2016-01-27 14:57:28Z spm $

N=length(x);

try, pmax=p; catch, pmax=min(ceil(N/10),256); end

for i=1:pmax,
    y(:,i)=x(pmax-i+1:end-i);
end
C=cov(y);
c=C(1,:);
gamma=c(2:end);
gamma0=c(1);
r=gamma/gamma0;

G=[];
for j=1:floor(pmax/2)-1,
    % Sum of adjacent pairs of autocovariances
    G(j)=gamma(2*j)+gamma(2*j+1);
end

if ~isempty(G)
    % Find minimum j such that all G's up to j are positive
    Gneg=find(G<0);
    if isempty(Gneg)
        m=length(G);
    else
        m1=min(Gneg);
        m=m1-1;
    end
else
    m=0;
end

ess=N/(1+2*sum(r(1:2*m)));

% figure;
% plot(c);
% title('Autocovariance');
% 
% figure
% plot(G);
% title('Sum of adjacent covariances');