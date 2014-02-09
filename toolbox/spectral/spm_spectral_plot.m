function spm_spectral_plot(Hz,csd,str,xlab,ylab)
% subplot for spectral arrays
% FORMAT spm_spectral_plot(Hz,csd,str,xlab,ylab)
%
% str  - format (default: '-')
% xlab - xlabel (default: 'Hz')
% ylab - ylabel (default: 'power')
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_spectral_plot.m 5873 2014-02-09 14:40:52Z karl $


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

