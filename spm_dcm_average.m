function spm_dcm_average (P,name)
% Produce an aggregate DCM model using Bayesian FFX averaging
% FORMAT spm_dcm_average (P,name)
%
% P         -  character/cell array of DCM filenames
% name      -  name of DCM output file (will be prefixed by 'DCM_avg_')
%
% This routine creates a new DCM model in which the parameters are averaged
% over a number of fitted DCM models. These can be over sessions or over
% subjects. This average model can then be interrogated using the standard
% DCM 'review' options to look at contrasts of parameters. The resulting
% inferences correspond to a Bayesian Fixed Effects analysis.
%
% Note that the Bayesian averaging is only applied to the A, B and C
% matrices (and matrix D if a nonlinear model is used).
% All other quantities in the average model are initially simply copied from
% the first DCM in the list. Subsequently, they are deleted before saving
% the average DCM in order to avoid any false impression that averaged
% models could be used for model comparison or contained averaged time series.
% Neither operation is valid and will be prevented by the DCM interface.
% Finally, note that only models with exactly the same A,B,C(,D) structure
% and the same brain regions can be averaged.
%
% A Bayesian random effects analysis can be implemented for a particular
% contrast using the spm_dcm_sessions.m function.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny & Klaas Enno Stephan
% $Id: spm_dcm_average.m 4185 2011-02-01 18:46:18Z guillaume $

try
    P;
catch
    [P, sts] = spm_select([2 Inf],'^DCM.*\.mat$','Select DCM*.mat files');
    if ~sts, return; end
end

try
    name;
catch
    name     = spm_input('Name for DCM_avg_???.mat','+1','s');
end

if ischar(P), P = cellstr(P); end
N = numel(P);

%-Loop through all selected models and get posterior means and precisions
%==========================================================================
for model = 1:N

    load(P{model});

    % Only look at those parameters with non-zero prior variance
    %----------------------------------------------------------------------
    pCdiag = diag(DCM.M.pC);
    wsel   = find(pCdiag);

    if model == 1
        wsel_first = wsel;
        DCM_first  = DCM;
    else
        if length(wsel) ~= length(wsel_first) || any(wsel ~= wsel_first)
            error('DCMs must have same structure.');
        end
    end

    % Get posterior precision matrix and mean
    %-------------------------------------------------------------------
    Cp              = DCM.Cp;
    Ep              = spm_vec(DCM.Ep);
    miCp(:,:,model) = inv(full(Cp(wsel,wsel)));
    mEp(:,model)    = full(Ep(wsel));

end


%-Average models using Bayesian fixed-effects analysis -> average Ep,Cp
%==========================================================================

% averaged posterior covariance
%--------------------------------------------------------------------------
Cp(wsel,wsel) = inv(sum(miCp,3));

% averaged posterior mean
%--------------------------------------------------------------------------
pE            = DCM.Ep;
wEp           = 0;
for model=1:N
    wEp       = wEp + miCp(:,:,model) * mEp(:,model);
end
Ep(wsel)      = Cp(wsel,wsel) * wEp;
Ep            = spm_unvec(Ep,pE);


%-Copy contents of first DCM into the output DCM and add BPA
%==========================================================================
DCM           = DCM_first;
DCM.models    = char(P);
DCM.averaged  = true;

% compute posterior probabilities and variance
%--------------------------------------------------------------------------
sw      = warning('off','SPM:negativeVariance');
Vp      = diag(Cp);
Pp      = 1 - spm_Ncdf(0,abs(spm_vec(Ep) - spm_vec(pE)),Vp);
warning(sw);

DCM.Ep  = Ep;
DCM.Cp  = Cp;
DCM.Vp  = spm_unvec(Vp,pE);
DCM.Pp  = spm_unvec(Pp,pE);


%-Save new DCM
%==========================================================================
DCM.name = [name ' (Bayesian FFX average)'];
if spm_check_version('matlab','7') >= 0
    save(['DCM_avg_' name '.mat'], 'DCM', '-V6');
else
    save(['DCM_avg_' name '.mat'], 'DCM');
end

% Warn the user how this average DCM should NOT be used
%--------------------------------------------------------------------------
disp(['Results of averaging DCMs were saved in DCM_avg_' name '.mat.']);
disp('Please note that this file only contains average parameter estimates');
disp('and their posterior probabilities, but NOT averaged time series.');
disp('Also, note that this file can NOT be used for model comparisons.');
