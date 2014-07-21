function spm_dcm_tfm_response(xY,pst,hz,top)
% displays evoked and induced responses
% FORMAT spm_dcm_tfm_response(xY,pst,hz)
%
% xY.erp{i} - (t x n):         an array over t time bins and n channels for
%                              condition i
% xY.csd{i} - (t x w x n x n): an array over t time bins, w frequency bins
%                              and n times n channels
%    pst - peristimulus time (seconds)
%    Hz  - frequency range   (Hz)
%__________________________________________________________________________
%
% This routine displays complex evoked and induced responses over peri-
% stimulus time in terms of 90% confidence intervals about the ERP and as 
% images of the spectral density for each cannel:
%
% see also spm_dcm_tfm_image - for between channel (coherence) responses)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_tfm_response.m 6112 2014-07-21 09:39:53Z karl $
 
% setup and defaults
%--------------------------------------------------------------------------
if nargin < 2, pst = 1:size(xY.csd{1},1); end
if nargin < 3, hz  = 1:size(xY.csd{1},2); end
if nargin < 4, top = 1;                   end
 
 
% plot time frequency responses
%==========================================================================
pst   = pst(:)';                                   % pst in ms
ne    = length(xY.csd);                            % number of event types
nc    = size(xY.erp{1},2);                         % number of channels
bands = kron([8; 13; 32],[1 1]);
for i = 1:nc
    for e = 1:ne
        j = 2*(i - 1)*ne + 2*(e - 1) + top;
        if j < 8
            % evoked response
            %--------------------------------------------------------------
            subplot(4,2,j)
            
            erp = xY.erp{e}(:,i)';
            csd = xY.csd{e}(:,:,i,i)';
           
            
            spm_plot_ci(erp,mean(abs(csd)),pst)
            str = sprintf('evoked: channel/source %i',i);
            title(str,'FontSize',16)
            xlabel('peristimulus time (ms)')
            spm_axis tight
            
            % induced response
            %--------------------------------------------------------------
            subplot(4,2,j + 1)
            csd = spm_detrend(abs(csd)')';
   
            imagesc(pst,hz,csd);
            str = sprintf('(adj.) induced: condition %i',e);
            title(str,'FontSize',16)
            xlabel('peristimulus time (ms)')
            ylabel('Hz'), axis xy
            
            % frequency ranges
            %--------------------------------------------------------------
            hold on; plot([pst(1) pst(end)],bands,':w'), hold off
            
        end
    end
end
drawnow
