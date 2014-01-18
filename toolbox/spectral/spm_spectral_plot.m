function spm_spectral_plot(Hz,csd,str,xlab,ylab)
% subplot for spectral arrays
% FORMAT spm_spectral_plot(Hz,csd,str,xlab,ylab)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_spectral_plot.m 5837 2014-01-18 18:38:07Z karl $


% order
%==========================================================================
m     = size(csd,2);

% plot
%--------------------------------------------------------------------------
for i = 1:m
    for j = 1:m
        subplot(m,m,(i - 1)*m + j)
        plot(Hz,abs(csd(:,i,j)),str), hold on
        xlabel(xlab)
        ylabel(ylab)
        axis square
    end
end

