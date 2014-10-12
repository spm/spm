function spm_dcm_tfm_response(xY,pst,hz,XY)
% displays evoked and induced responses
% FORMAT spm_dcm_tfm_response(xY,pst,hz,[XY])
%
% xY.erp{i} - (t x n):         an array over t time bins and n channels for
%                              condition i
% xY.csd{i} - (t x w x n x n): an array over t time bins, w frequency bins
%                              and n times n channels
%    pst - peristimulus time (seconds)
%    Hz  - frequency range   (Hz)
%
% XY true value for overplotting
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
% $Id: spm_dcm_tfm_response.m 6234 2014-10-12 09:59:10Z karl $

% setup and defaults
%--------------------------------------------------------------------------
if nargin < 2, pst = 1:size(xY.erp{1},1); end
if nargin < 3, hz  = 1:size(xY.csd{1},2); end


% plot time frequency responses
%==========================================================================
pst   = pst(:)';                                   % pst in ms
ne    = length(xY.csd);                            % number of event types
nc    = size(xY.erp{1},2);                         % number of channels
bands = kron([8; 13; 32],[1 1]);

if nargin < 4
    for i = 1:nc
        for e = 1:ne
            j = 2*(i - 1)*ne + 2*(e - 1) + 1;
            if j < 8
                
                % evoked response
                %----------------------------------------------------------
                subplot(4,2,j)
                
                erp = xY.erp{e}(:,i)';
                csd = xY.csd{e}(:,:,i,i)';
                
                
                spm_plot_ci(erp,mean(abs(csd)),pst)
                str = sprintf('evoked: channel/source %i',i);
                title(str,'FontSize',16)
                xlabel('peristimulus time (ms)')
                spm_axis tight
                
                % induced response
                %----------------------------------------------------------
                subplot(4,2,j + 1)
                csd = spm_detrend(abs(csd)')';
                
                imagesc(pst,hz,csd);
                str = sprintf('(adj.) induced: condition %i',e);
                title(str,'FontSize',16)
                xlabel('peristimulus time (ms)')
                ylabel('Hz'), axis xy
                
                % frequency ranges
                %----------------------------------------------------------
                hold on; plot([pst(1) pst(end)],bands,':w'), hold off
                
            end
        end
    end
    
else
    
    for i = 1:nc
        for e = 1:ne
            j = 3*(i - 1)*ne + 3*(e - 1) + 1;
            if j < 12
                
                % evoked response
                %----------------------------------------------------------
                subplot(4,3,j)
                
                erp = xY.erp{e}(:,i)';
                csd = xY.csd{e}(:,:,i,i)';
                ERP = XY.erp{e}(:,i)';
                CSD = XY.csd{e}(:,:,i,i)';
                
                
                spm_plot_ci(erp,mean(abs(csd)),pst), hold on
                plot(pst,ERP,':'), hold off
                str = sprintf('evoked: channel/source %i',i);
                title(str,'FontSize',16)
                xlabel('peristimulus time (ms)')
                spm_axis tight
                
                % induced response
                %----------------------------------------------------------
                subplot(4,3,j + 1)
                csd = spm_detrend(abs(csd)')';
                
                imagesc(pst,hz,csd);
                str = sprintf('(adj.) induced: condition %i',e);
                title(str,'FontSize',16)
                xlabel('peristimulus time (ms)')
                ylabel('Hz'), axis xy
                
                % frequency ranges
                %----------------------------------------------------------
                hold on; plot([pst(1) pst(end)],bands,':w'), hold off
                
                % induced response (observed)
                %----------------------------------------------------------
                subplot(4,3,j + 2)
                CSD = spm_detrend(abs(CSD)')';
                
                imagesc(pst,hz,CSD);
                str = sprintf(' - empirical',e);
                title(str,'FontSize',16)
                xlabel('peristimulus time (ms)')
                ylabel('Hz'), axis xy
                
                % frequency ranges
                %----------------------------------------------------------
                hold on; plot([pst(1) pst(end)],bands,':w'), hold off
                
            end
        end
    end
    drawnow
end
