function spm_spectral_plot(Hz,csd,str,xlab,ylab)
% subplot for spectral arrays
% FORMAT spm_spectral_plot(Hz,csd,str,xlab,ylab)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_spectral_plot.m 5853 2014-01-24 20:38:11Z karl $


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
        if i == j
            title('auto','FontSize',16)
        elseif j > i
            title('backward','FontSize',16)
        elseif j < i
            title('forward','FontSize',16)
        end
        axis square
    end
end

