function [] = mci_compare_sensitivities (model,pars)
% Compare methods for sensitivity computation
% FORMAT [] = mci_compare_sensitivities (model,pars)
%
% model     'phase', 'nmm-r2p2'
% pars      vector indicating which sensitivities to plot
%           eg. [1,2,..,Np] (default) for all parameters
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: mci_compare_sensitivities.m 6548 2015-09-11 12:39:47Z will $

[P,M,U,Y] = mci_compare_setup (model);

if nargin < 2 | isempty(pars)
    pars=[1:M.Np];
end
Np=length(pars);

disp(' ');
disp('Plot sensitivity to parameter changes as computed by:');
disp('(1) Forward equations based on Matlab''s ODE suite (red)');
disp('(2) Forward equations based on Sundials (blue)');

[G,sy_m,st] = spm_mci_sens (P,M,U);

[G,sy,st] = spm_mci_sens_sun (P,M,U);

if st==-1
    disp('Problem with integration');
    return
end

mci_plot_outputs(M,G);

hs=figure;
set(hs,'Name','Sensitivities');
k=1; lw=2;
for i=1:M.l,
    for p=1:Np,
        j=pars(p);
        subplot(M.l,Np,k);
        plot(M.t,squeeze(sy(:,i,j)),'LineWidth',lw);
        hold on
        plot(M.t,squeeze(sy_m(:,i,j)),'r','LineWidth',lw);
        grid on
        title (sprintf('dy(%d)/dp(%d)',i,j));
        xlabel('Time');
        k=k+1;
    end
end
