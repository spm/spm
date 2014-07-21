function spm_dcm_tfm_transfer(dtf,pst,Hz)
% displays time-frequency modulation transfer functions
% FORMAT spm_dcm_tfm_transfer(dtf,pst,Hz)
% 
% dtf - (t x w x n x u): an array over t time bins, w frequency bins,
%                        n channels and u inputs
% pst - peristimulus time (seconds)
% Hz  - frequency range (Hz)
%__________________________________________________________________________
%
% This routine displays complex modulation transfer functions over 
% peristimulus time as images of the absolute values and first order
% kernels mapping from endogenous inputs to neuronal states
%
% See also: spm_dcm_tfm_responses (and spm_dcm_tfm_image) 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_tfm_transfer.m 6112 2014-07-21 09:39:53Z karl $
 
% setup and defaults
%--------------------------------------------------------------------------
if nargin < 2, pst = 1:size(dtf,1); end
if nargin < 3, Hz  = 1:size(dtf,2); end


% plot modulation transfer functions (in frequency space) - MTF
%==========================================================================
nc    = size(dtf,3);
nu    = size(dtf,4);
bands = kron([8; 13; 32],[1 1]);
for i = 1:nc
    for j = 1:nu
        
        subplot(2*nc,nu,(i - 1)*nc + j)
        imagesc(pst,Hz,abs(dtf(:,:,i,j)).^2');
        str  = sprintf('transfer function: %i to %i',j,i);
        title(str)
        xlabel('peristimulus time (ms)')
        ylabel('Hz'), axis xy
        hold on; plot([pst(1) pst(end)],bands,':w'), hold off
        
    end
end
 
% plot transfer functions (impulse response function) - IRF
%==========================================================================

% get inverse Fourier transform of MTF
%--------------------------------------------------------------------------
[irf, lag]  = spm_csd2ccf(dtf,Hz);
lag   = lag*1000;

% restrict range of lags
%--------------------------------------------------------------------------
j     = find(0 < lag & lag < 256);
lag   = lag(j);
irf   = irf(:,j,:,:);

for i = 1:nc
    for j = 1:nu
        
        subplot(2*nc,nu, nc*nu + (i - 1)*nc + j)
        plot(lag,real(irf(:,:,i,j)));
        str  = sprintf('kernel: %i to %i',j,i);
        title(str)
        xlabel('lag (ms)')
        spm_axis tight
        
    end
end
drawnow
