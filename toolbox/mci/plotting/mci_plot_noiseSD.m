function [] = mci_plot_noiseSD (Ce,ind)
% Plot posterior over observation noise SDs
% FORMAT [] = mci_plot_noiseSD (Ce,ind)
%
% Ce    [d x d x Ns] d-dimensional observations with Ns samples
% ind   indices of samples to plot
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: mci_plot_noiseSD.m 6548 2015-09-11 12:39:47Z will $

d=size(Ce,1);
dist.ind=ind;
for j=1:d,
    dist.P(j,:)=sqrt(squeeze(Ce(j,j,:))');
    dist.names{j}=sprintf('NoiseSD (%d)',j);
end
dist.color='k';
dist.type='sample';
dist.ks=1;
figure
rd=ceil(sqrt(d));
for j=1:d,
    subplot(rd,rd,j);
    mci_plot_dist(dist,j);
    grid on
end
