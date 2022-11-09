function [pC] = spm_dcm_lock(pV)
% Lock experimental effects by introducing prior correlations
% FORMAT [pC] = spm_dcm_lock(pV)
%__________________________________________________________________________
%
% pV   - prior variance
% pC   - prior covariance
%__________________________________________________________________________

 
% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging

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
