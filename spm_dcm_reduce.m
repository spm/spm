function DCM = spm_dcm_reduce(DCM,rE,rC)
% Reduce the posterior of DCM given new priors (rE,rC)
% FORMAT RCM = spm_dcm_reduce(DCM,rE,rC)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2015-2022 Wellcome Centre for Human Neuroimaging


% deal with cell arrays
%--------------------------------------------------------------------------
if iscell(DCM)
    for i = 1:numel(DCM)
       DCM{i} = spm_dcm_reduce(DCM{i},rE,rC);
    end
    return
end

% empirical prior and posterior densities
%--------------------------------------------------------------------------
pE = DCM.M.pE;
pC = DCM.M.pC;
qE = DCM.Ep;
qC = DCM.Cp;

% evaluate posteriors under original priors
%--------------------------------------------------------------------------
[F,sE,sC] = spm_log_evidence_reduce(qE,qC,pE,pC,rE,rC);

DCM.M.pE = rE;
DCM.M.pC = rC;
DCM.Ep   = sE;
DCM.Cp   = sC;
DCM.F    = F + DCM.F;
