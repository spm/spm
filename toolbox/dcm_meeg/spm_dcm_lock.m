function [pC] = spm_dcm_lock(pV)
% locks experimental effects by introducing prior correlations
% FORMAT [pC] = spm_dcm_lock(pV)
%__________________________________________________________________________
%
% pV   - prior variance
% pC   - prior covariance
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_lock.m 4095 2010-10-22 19:37:51Z karl $

% lock experimental effects by introducing prior correlations
%==========================================================================
pC    = spm_diag(spm_vec(pV));
for i = 1:length(pV.B)
    pB      = pV;
    pB.B{i} = pB.B{i} - pB.B{i};
    pB      = spm_vec(pV)  - spm_vec(pB);
    pB      = sqrt(pB*pB') - diag(pB);
    pC      = pC + pB;
end
