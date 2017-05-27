function [BPA] = spm_dcm_bpa(P,nocd)
% Produce an aggregate DCM using Bayesian parameter averaging
% FORMAT [BPA] = spm_dcm_bpa(DCM,nocd)
%
% DCM  - {N [x M]} structure array of DCMs from N subjects
% ------------------------------------------------------------
%     DCM{i}.M.pE - prior expectations of P parameters
%     DCM{i}.M.pC - prior covariance
%     DCM{i}.Ep   - posterior expectations
%     DCM{i}.Cp   - posterior covariance
%
% nocd - optional flag for suppressing conditional dependencies.
%        This is useful when evaluating the BPA of individual (contrasts
%        of) parameters, where the BPA of a contrast should not be confused
%        with the contrast of a BPA.
%
% BPA  - DCM structure (array) containing Bayesian parameter averages
% ------------------------------------------------------------
%     BPA.M.pE - prior expectations of P parameters
%     BPA.M.pC - prior covariance
%     BPA.Ep   - posterior expectations
%     BPA.Cp   - posterior covariance
%
%     BPA.Pp   - posterior probability of > 0
%     BPA.Vp   - posterior variance
%     BPA....  - other fields from DCM{1[,:]}
%__________________________________________________________________________
%
% This routine creates a new DCM in which the parameters are averaged over
% a number of fitted DCMs. These can be over sessions or over subjects.
% This average model can then be interrogated using the standard DCM
% 'review' options to look at contrasts of parameters. The resulting
% inferences correspond to a Bayesian Fixed Effects analysis. If called
% with no output arguments the Bayesian parameter average DCM will be
% written to DCM_BPA.mat; otherwise, the DCM structure is returned as BPA.
% 
% If DCM is an {N x M} array, Bayesian parameter averaging will be
% applied to each model (i.e., each row) - and BPA becomes a {1 x M} cell 
% array.
%
% See also spm_dcm_bma.m, spm_dcm_bmr.m and spm_dcm_peb.m
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny & Klaas Enno Stephan
% $Id: spm_dcm_bpa.m 7081 2017-05-27 19:36:09Z karl $

 
% Preiminaries
%--------------------------------------------------------------------------
if nargin < 1
    [P, sts] = spm_select([1 Inf],'^DCM.*\.mat$','Select DCM*.mat files');
    if ~sts, return; end
end
if ischar(P), P = cellstr(P); end

% conditional dependencies option
%--------------------------------------------------------------------------
if nargin < 2, nocd = 0; end

% repeat for each column if P is and array
%--------------------------------------------------------------------------
if ~isvector(P)
    
    % loop over models in each column
    %----------------------------------------------------------------------
    for i = 1:size(P,2)
        BPA{i}  = spm_dcm_bpa(P(:,i),nocd);
    end
    
    return
end

 
%-Loop through all selected models and get posterior means and precisions
%==========================================================================
N     = numel(P);
TOL   = exp(-16);
sEp   = []; 
Up    = []; 
Pp    = [];
pE    = [];

for i = 1:N
    
    % get DCM structure
    %----------------------------------------------------------------------
    if ischar(P{i}), load(P{i}); else DCM = P{i}; end
    
    % Only look at those parameters with non-zero prior variance
    %----------------------------------------------------------------------
    if isstruct(DCM.M.pC)
        pC = diag(spm_vec(DCM.M.pC));
    else
        pC = DCM.M.pC;
    end
    j   = diag(pC) > TOL & diag(pC) < 1/TOL;
    pC  = diag(j)*pC*diag(j);
    
    % find the space spanned by the prior covariance
    %----------------------------------------------------------------------
    if i == 1
        
        % suppress prior dependencies if necessary
        %----------------------------------------------------------------------
        if nocd, pC = diag(diag(pC)); end
        
        U       = spm_svd(pC,0);
        j_first = j;
        BPA     = DCM;
        
    else
        if any(j ~= j_first)
            error('DCMs must have same structure.');
        end
    end
    
    % suppress posterior dependencies if necessary
    %----------------------------------------------------------------------
    Ep          = spm_vec(DCM.Ep);
    Cp          = DCM.Cp;
    if nocd, Cp = diag(diag(Cp)); end
    
    % Get posterior precision matrix and mean
    %----------------------------------------------------------------------
    sEp(:,i)    = Ep;    
    Up(:,i)     = U'*Ep;
    Cp          = U'*Cp*U;
    Pp(:,:,i)   = inv(full(Cp));
 
end

 
%-Average models using Bayesian fixed-effects analysis -> average Ep,Cp
%==========================================================================

% average posterior covariance and mean
%--------------------------------------------------------------------------
Cp  = inv(sum(Pp,3));
mEp = mean(sEp,2);

% average posterior mean
%--------------------------------------------------------------------------
Ep    = 0;
for i = 1:N
    Ep = Ep + Pp(:,:,i)*Up(:,i);
end
Ep  = Cp*Ep;

% project back through U
%--------------------------------------------------------------------------
Cp  = U*Cp*U';
Ep  = U*Ep + mEp - U*(U'*mEp);
Ep  = spm_unvec(Ep,DCM.M.pE);

 
%-Copy contents of first DCM into the output DCM and add BPA
%==========================================================================
BPA.averaged   = true;
try
    BPA.models = char(P);
end

% compute posterior probabilities and variance
%--------------------------------------------------------------------------
sw     = warning('off','SPM:negativeVariance');
Vp     = diag(Cp);
alpha  = 0;
Pp     = 1 - spm_Ncdf(0,abs(spm_vec(Ep) - alpha),Vp);
warning(sw);

BPA.Ep = Ep;
BPA.Cp = Cp;
BPA.Vp = spm_unvec(Vp,Ep);
BPA.Pp = spm_unvec(Pp,Ep);

 
%-Save new DCM if there are no output arguments
%==========================================================================
if nocd
    BPA.name = 'BPA - no conditional dependencies';
else
    BPA.name = 'Bayesian parameter average';
end
if ~nargout
    save('DCM_BPA.mat', 'BPA', spm_get_defaults('mat.format'));
end