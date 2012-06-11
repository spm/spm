function spm_dcm_tfm_image(csd,top,pst,hz)
% displays time-frequency complex cross spectra
% FORMAT spm_dcm_tfm_image(csd,top,pst,hz)
% 
% csd - (t x w x n x n): a data array over t time bins, w frequency bins
%                       and n times n channels
% top - switch to display at the top or bottom of the current figure
% pst - peristimulus time (for plotting)
% Hz  - frequency range (for plotting)
%__________________________________________________________________________
%
% this routine displays complex cross spectra over peristimulus time as
% images of the absolute values (coherence) and cross covariance functions
% over pairs of channels
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_tfm_image.m 4768 2012-06-11 17:06:55Z karl $
 
% setup and defaults
%--------------------------------------------------------------------------
if nargin < 2, top = 1;           end
if nargin < 3, pst = 1:size(csd,1); end
if nargin < 4, hz  = 1:size(csd,2); end


% plot time frequency responses
%==========================================================================
nc    = size(csd,3);
bands = kron([8; 13; 32],[1 1]);
for i = 1:nc
    
    % evaluate cross covariance function
    %----------------------------------------------------------------------
    g         = abs(csd);
    [ccf,lag] = spm_csd2ccf(csd,hz);
    lag       = lag*1000;
    j         = find(-64 < lag & lag < 64);
    lag       = lag(j);
    ccf       = ccf(:,j,:,:);
    
    
    % spectral power
    %----------------------------------------------------------------------
    for j = i:i
        subplot(2*nc,nc,~top*nc*nc + (i - 1)*nc + j)
        
        imagesc(pst,hz,abs(g(:,:,i,j)).^2');
        title('Spectral density','FontSize',16)
        xlabel('peristimulus time (ms)')
        ylabel('Hz'), axis xy
        
        hold on; plot([pst(1) pst(end)],bands,':w'), hold off
        
    end
    
    % coherence functions
    %----------------------------------------------------------------------
    for j = (i + 1):nc
        subplot(2*nc,nc,~top*nc*nc + (i - 1)*nc + j)
        
        imagesc(pst,hz,abs(g(:,:,i,j)).^2');
        title('Coherence','FontSize',16)
        xlabel('peristimulus time (ms)')
        ylabel('Hz'), axis xy
        
        hold on; plot([pst(1) pst(end)],bands,':w'), hold off
        
    end
    
    % cross covariance functions
    %----------------------------------------------------------------------
    for j = 1:(i - 1)
        subplot(2*nc,nc,~top*nc*nc + (i - 1)*nc + j)
        
        imagesc(pst,lag,ccf(:,:,i,j)');
        title('Cross-covariance','FontSize',16)
        xlabel('peristimulus (ms)')
        ylabel('lag (ms)'), axis xy
        
        hold on; plot([pst(1) pst(end)],[0 0],'-.w'), hold off
        
    end
end
drawnow
