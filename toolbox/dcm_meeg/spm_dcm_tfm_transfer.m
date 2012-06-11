function spm_dcm_tfm_transfer(mtf,pst,hz)
% displays time-frequency transfer functions (in time and frequency space)
% FORMAT spm_dcm_tfm_transfer(mtf,pst,hz)
% 
% mtf - (t x w x n x u): a MTF array over t time bins, w frequency bins,
%                        n channels and u inputs
% pst - peristimulus time (for plotting)
% Hz  - frequency range (for plotting)
%__________________________________________________________________________
%
% This routine displays complex modulation transfer functions over 
% peristimulus time as images of the absolute values and first order
% kernels mapping from endogenous inputs to neuronal states
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_tfm_transfer.m 4768 2012-06-11 17:06:55Z karl $
 
% setup and defaults
%--------------------------------------------------------------------------
if nargin < 2, pst = 1:size(mtf,1); end
if nargin < 3, hz  = 1:size(mtf,2); end
 
 
% plot modulation transfer functions (in frequency space) - MTF
%==========================================================================
nc    = size(mtf,3);
nu    = size(mtf,4);
bands = kron([8; 13; 32],[1 1]);
for i = 1:nc
    for j = 1:nu
        subplot(2*nc,nu,(i - 1)*nc + j)
        
        imagesc(pst,hz,abs(mtf(:,:,i,j)).^2');
        str  = sprintf('tansfer function: %i to %i',j,i);
        title(str,'FontSize',16)
        xlabel('peristimulus time (ms)')
        ylabel('Hz'), axis xy
        
        hold on; plot([pst(1) pst(end)],bands,':w'), hold off
        
    end
end
 
% plot transfer functions (impulse response function) - IRF
%==========================================================================

% get inverse Fourier transform of MTF
%--------------------------------------------------------------------------
[irf, lag]  = spm_csd2ccf(mtf,hz);

% restrict range of lags
%--------------------------------------------------------------------------
lag   = lag*1000;
j     = find(0 < lag & lag < 128);
lag   = lag(j);
irf   = irf(:,j,:,:);

for i = 1:nc
    for j = 1:nu
        subplot(2*nc,nu, nc*nu + (i - 1)*nc + j)
        
        plot(lag,real(irf(:,:,i,j)));
        str  = sprintf('kernel: %i to %i',j,i);
        title(str,'FontSize',16)
        xlabel('lag (ms)')
        
    end
end
drawnow
