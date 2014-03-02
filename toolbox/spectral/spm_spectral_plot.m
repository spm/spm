function spm_spectral_plot(Hz,csd,str,xlab,ylab,r)
% subplot for spectral arrays
% FORMAT spm_spectral_plot(Hz,csd,str,xlab,ylab)
%
% str  - format (default: '-')
% xlab - xlabel (default: 'Hz')
% ylab - ylabel (default: 'power')
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_spectral_plot.m 5902 2014-03-02 18:26:06Z karl $


% order
%==========================================================================
m     = size(csd,2);

if nargin < 3, str  = '-'; end
if nargin < 4, xlab = 'Frequency'; end
if nargin < 5, ylab = 'Power'; end

% plot
%--------------------------------------------------------------------------
for i = 1:m
    for j = 1:m
        subplot(m,m,(i - 1)*m + j)
        if nargin > 5
            plot(Hz,real(csd(:,i,j)),str), hold on
            plot(Hz,imag(csd(:,i,j)),':'), hold on
        else
            plot(Hz,abs(csd(:,i,j)),str), hold on
        end
        xlabel(xlab)
        ylabel(ylab)
        if i == j
            title('auto','FontSize',16)
        elseif j > i
            title('backward','FontSize',16)
        elseif j < i
            title('forward','FontSize',16)
        end
        axis square, set(gca,'XLim',[min(Hz(1),0) Hz(end)])
    end
end

