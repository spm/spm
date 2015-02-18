function [i,pC,pE] = spm_find_pC(pC,pE,fields)
% Utility routine that finds the indices of non-zero covariance
% FORMAT [i]       = spm_find_pC(pC,pE,fields])
% FORMAT [i,rC,rE] = spm_find_pC(DCM)
% 
% pC     - covaraince matrix or variance stucture
% pE     - parameter structure
% fields - desired feilds of pE
%
% or
%
% DCM    - DCM structure
%
% i      - find(diag(pC) > TOL)
% rC     - reduced covariances
% rE     - reduced expectation
% 
%__________________________________________________________________________
% Copyright (C) 2006-2013 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_find_pC.m 6342 2015-02-18 15:22:07Z karl $


%-get rC
%--------------------------------------------------------------------------
if isfield(pC,'options')
    [pC,pE] = spm_find_rC(pC);
end

%-Deal with structures
%--------------------------------------------------------------------------
if isstruct(pC)
    q = spm_vec(pC);
else
    q = diag(pC);
end

%-Get indices
%--------------------------------------------------------------------------
i  = find(q > mean(q(q < 1024))/1024);

%-subsample fields if necessary
%--------------------------------------------------------------------------
if nargin > 2
    if isstruct(pE)
        j = spm_fieldindices(pE,fields{:});
        if ~isempty(j)
            i = j(ismember(j,i));
        end
    end
end

return


function [pC,pE] = spm_find_rC(DCM)
% FORMAT:[pC,pE] = spm_find_rC(DCM)
% model priors
%__________________________________________________________________________

% Get full priors and posteriors
%--------------------------------------------------------------------------
if isfield(DCM,'M')
    pC  = DCM.M.pC;
    pE  = DCM.M.pE;
    return
end

% get priors from model specification
%------------------------------------------------------------------
if isfield(options,'analysis')
    if strcmpi(options.analysis,'IND')
        [pE,~,pC] = spm_ind_priors(DCM.A,DCM.B,DCM.C,Nf);
    else
        [pE,pC] = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,options.model);
    end
else
    [pE,pC] = spm_dcm_fmri_priors(DCM.a,DCM.b,DCM.c,DCM.d,options);
end
