function [coh,fsd] = spm_csd2coh(csd,Hz)
% Converts cross spectral density to coherence and (phase) delay
% FORMAT [coh,fsd] = spm_csd2coh(csd,Hz)
%
% csd  (Hz,:,:) - cross spectral density (cf, mar.P)
% Hz   (n x 1)  - vector of frequencies
%
% coh           - coherence
% fsd           - frequency specific delay (seconds) 
%               - phase-delay/radial frequnecy
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_csd2coh.m 5816 2013-12-23 18:52:56Z karl $
 

% unpack cells
%--------------------------------------------------------------------------
if iscell(csd)
    for i = 1:length(csd)
       [cohi,fsdi] = spm_csd2coh(csd{i},Hz);
       coh{i}     = cohi;
       fsd{i}     = fsdi;
    end
    return
end

% preliminaries
%--------------------------------------------------------------------------
Hz    = spm_vec(Hz);
 
% compute coherence and delay
%==========================================================================
for i = 1:size(csd,2)
    for j = 1:size(csd,2)
        coh(:,i,j) = abs(csd(:,i,j).*conj(csd(:,i,j)))./abs(csd(:,i,i).*csd(:,j,j));
        fsd(:,i,j) = unwrap(angle(csd(:,i,j)))./(2*pi*Hz);
    end
end
