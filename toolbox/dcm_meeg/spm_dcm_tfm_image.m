function spm_dcm_tfm_image(csd,pst,hz,top)
% displays time-frequency complex cross spectra
% FORMAT spm_dcm_tfm_image(csd,top,pst,hz)
% 
% csd - (t x w x n x n): a data array over t time bins, w frequency bins
%                       and n times n channels
% pst - peristimulus time (for plotting)
% Hz  - frequency range (for plotting)
% top - [0/1] switch to display at the top or bottom of the current figure
%__________________________________________________________________________
%
% This routine displays complex cross spectra over peristimulus time as
% images of the absolute values (coherence) and cross covariance functions
% over pairs of channels.
%
% See also: spm_dcm_tfm_responses (and spm_dcm_tfm_transfer) 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_tfm_image.m 6112 2014-07-21 09:39:53Z karl $
 
% setup and defaults
%--------------------------------------------------------------------------
two = 2;
if nargin < 2, pst = 1:size(csd,1); end
if nargin < 3, hz  = 1:size(csd,2); end
if nargin < 4, top = 0; two = 1;    end


% plot time frequency responses
%==========================================================================
nc    = size(csd,3);
bands = kron([8; 13; 32],[1 1]);


% evaluate cross covariance function
%--------------------------------------------------------------------------
[ccf,lag] = spm_csd2ccf(csd,hz,1/1024);
lag       = lag*1000;
j         = find(-64 < lag & lag < 64);
lag       = lag(j);
ccf       = ccf(:,j,:,:);


for i = 1:nc    
    
    % spectral power
    %----------------------------------------------------------------------
    for j = i:i
        
        subplot(two*nc,nc,top*nc*nc + (i - 1)*nc + j)
        imagesc(pst,hz,abs(csd(:,:,i,j)).^2');
        title('Spectral density')
        xlabel('peristimulus time (ms)')
        ylabel('Hz'), axis xy
        hold on; plot([pst(1) pst(end)],bands,':w'), hold off
        
    end
    
    % coherence functions
    %----------------------------------------------------------------------
    for j = (i + 1):nc
        
        subplot(two*nc,nc,top*nc*nc + (i - 1)*nc + j)
        imagesc(pst,hz,abs(csd(:,:,i,j)).^2');
        title('Coherence')
        xlabel('peristimulus time (ms)')
        ylabel('Hz'), axis xy
        hold on; plot([pst(1) pst(end)],bands,':w'), hold off
        
    end
    
    % cross covariance functions
    %----------------------------------------------------------------------
    for j = 1:(i - 1)
        
        subplot(two*nc,nc,top*nc*nc + (i - 1)*nc + j)
        imagesc(pst,lag,ccf(:,:,i,j)');
        title('Cross-covariance')
        xlabel('peristimulus (ms)')
        ylabel('lag (ms)'), axis xy     
        hold on; plot([pst(1) pst(end)],[0 0],'-.k'), hold off
        
    end
end
drawnow
