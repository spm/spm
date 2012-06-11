function spm_dcm_tfm_response(xY,pst,hz)
% displays time-frequency complex cross spectra
% FORMAT spm_dcm_tfm_response(xY,pst,hz)
%
% xY.erp{i} - (t x n):         an array over t time bins and n channels for
%                              condition i
% xY.csd{i} - (t x w x n x n): an array over t time bins, w frequency bins
%                              and n times n channels
%    pst - peristimulus time (for plotting)
%    Hz  - frequency range   (for plotting)
%__________________________________________________________________________
%
% this routine displays complex cross spectra over peristimulus time as
% images of the absolute values (coherence) and cross covariance functions
% over pairs of channels
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_tfm_response.m 4768 2012-06-11 17:06:55Z karl $

% setup and defaults
%--------------------------------------------------------------------------
if nargin < 2, pst = 1:size(xY.csd{1},1); end
if nargin < 3, hz  = 1:size(xY.csd{1},2); end


% plot time frequency responses
%==========================================================================
ne    = length(xY.csd);                            % number of event types
nc    = size(xY.csd{1},3);                         % numberof channels
bands = kron([8; 13; 32],[1 1]);
for i = 1:nc
    for e = 1:ne
        try
            % evoked response
            %--------------------------------------------------------------
            subplot(4,2,2*(i - 1)*ne + 2*(e - 1) + 1)
            
            erp = xY.erp{e}(:,i)';
            csd = xY.csd{e}(:,:,i,i)';
            
            spm_plot_ci(pst,erp,mean(csd))
            str = sprintf('evoked: channel %i',i);
            title(str,'FontSize',16)
            xlabel('peristimulus time (ms)')
            spm_axis tight
            
            % induced response
            %--------------------------------------------------------------
            subplot(4,2,2*(i - 1)*ne + 2*(e - 1) + 2)
            
            imagesc(pst,hz,abs(csd).^2);
            str = sprintf('induced: condition %i',e);
            title(str,'FontSize',16)
            xlabel('peristimulus time (ms)')
            ylabel('Hz'), axis xy
            
            % frequnecy ranges
            %--------------------------------------------------------------
            hold on; plot([pst(1) pst(end)],bands,':w'), hold off
            
        end
    end
end
drawnow
