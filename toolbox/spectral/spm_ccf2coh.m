function [coh] = spm_ccf2coh(ccf,Hz)
% Converts cross covariance function to coherence
% FORMAT [coh] = spm_ccf2coh(ccf,Hz)
%
% ccf  (N,:,:)  - cross covariance functions
% Hz   (n x 1)  - vector of frequencies (Hz)
%
% coh           - coherence
%
% See also: spm_???2???.m
%     ??? = {'ccf','csd','gew','mar','coh','mtf','ker','ssm','dcm'}
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_ccf2coh.m 7774 2020-01-25 18:07:03Z karl $
 
% unpack cells
%--------------------------------------------------------------------------
if iscell(ccf)
    for i = 1:length(ccf)
       coh{i}    = spm_ccf2coh(ccf{i},Hz);
    end
    return
end

% convert via cross spectral density
%--------------------------------------------------------------------------
csd = spm_ccf2csd(ccf,Hz);
coh = spm_csd2coh(csd,Hz);
