function [DCM] = spm_dcm_average(P,name,nocond,graphics)
% Produce an aggregate DCM model using Bayesian FFX averaging
% FORMAT [DCM] = spm_dcm_average(P,name,nocond,graphics)
%
% P         -  character/cell array of DCM filenames
% name      -  name of DCM output file (will be prefixed by 'DCM_avg_')
% nocond    -  optional flag for suppressing conditional dependencies
% graphics  -  optional flag for showing outliers (based on conditional
%              entropy)
%
% This routine creates a new DCM in which the parameters are averaged
% over a number of fitted DCM models. These can be over sessions or over
% subjects. This average model can then be interrogated using the standard
% DCM 'review' options to look at contrasts of parameters. The resulting
% inferences correspond to a Bayesian Fixed Effects analysis. If called with
% no output arguments the Bayesian parameter average DCM will be written to
% file, otherwise the DCM structure is returned.
%
% Note that the Bayesian averaging is only applied to the A, B and C
% matrices (and matrix D if a nonlinear model is used).
% All other quantities in the average model are initially simply copied from
% the first DCM in the list. Subsequently, they are deleted before saving
% the average DCM in order to avoid any false impression that averaged
% models could be used for model comparison or contained averaged time series.
% Neither operation is valid and will be prevented by the DCM interface.
% Finally, note that only models with exactly the same A,B,C,(D) structure
% and the same brain regions can be averaged.
%
% A Bayesian random effects analysis can be implemented for a particular
% contrast using the spm_dcm_sessions.m function.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny & Klaas Enno Stephan
% $Id: spm_dcm_average.m 6724 2016-02-19 19:13:07Z karl $


% Preiminaries
%--------------------------------------------------------------------------
try
    P;
catch
    [P, sts] = spm_select([1 Inf],'^DCM.*\.mat$','Select DCM*.mat files');
    if ~sts, return; end
end
if ischar(P), P = cellstr(P); end
N = numel(P);

% options and filename
%--------------------------------------------------------------------------
try, name;     catch, name = spm_input('Name for DCM_avg_???.mat','+1','s'); end
try, nocond;   catch, nocond = 0; end
try, graphics; catch, graphics = 0; end


%-Loop through all selected models and get posterior means and precisions
%==========================================================================
for model = 1:N
    
    % get DCM structure
    %----------------------------------------------------------------------
    if ischar(P{model}), load(P{model}); else, DCM = P{model}; end
    
    % Only look at those parameters with non-zero prior variance
    %----------------------------------------------------------------------
    if isstruct(DCM.M.pC)
        pC = diag(spm_vec(DCM.M.pC));
    else
        pC = DCM.M.pC;
    end
    wsel   = diag(pC) > 0 & diag(pC) < 128;
    pC     = diag(wsel)*pC*diag(wsel);
    wsel   = find(wsel);
    
    % find the space spanned by the prior covariance
    %----------------------------------------------------------------------
    if model == 1
        
        % suppress prior dependencies if necessary
        %----------------------------------------------------------------------
        if nocond, pC = diag(diag(pC)); end
        
        U          = spm_svd(pC,0);
        wsel_first = wsel;
        DCM_first  = DCM;
    else
        if length(wsel) ~= length(wsel_first) || any(wsel ~= wsel_first)
            error('DCMs must have same structure.');
        end
    end
    
    % suppress posterior dependencies if necessary
    %----------------------------------------------------------------------
    Cp              = DCM.Cp;
    if nocond, Cp   = diag(diag(Cp)); end
    
    % Get posterior precision matrix and mean
    %----------------------------------------------------------------------
    Cp              = U'*Cp*U;
    Ep              = U'*spm_vec(DCM.Ep);
    miCp(:,:,model) = inv(full(Cp));
    mEp(:,model)    = Ep;
    
    
    % evaluate diagnostics
    %----------------------------------------------------------------------
    if graphics
        T(model) = trace(Cp);
        H(model) = spm_logdet(miCp(:,:,model));
        Fs(model) = DCM.F;
    end
    
end


%-Report free energies and conditional entropies
%==========================================================================
if graphics
    spm_figure('GetWin','BPA');
    
    subplot(3,1,1), bar(Fs)
    title('Free energy','FontSize',16)
    xlabel('Subject'), axis square
    
    subplot(3,1,2), bar(H)
    title('Conditional entropy','FontSize',16)
    xlabel('Subject'), axis square
    
    subplot(3,1,3), bar(T)
    title('Posterior variance','FontSize',16)
    xlabel('Subject'), axis square
end


%-Average models using Bayesian fixed-effects analysis -> average Ep,Cp
%==========================================================================

% averaged posterior covariance
%--------------------------------------------------------------------------
ipC = inv(U'*pC*U);
Cp  = inv(sum(miCp,3));

% averaged posterior mean
%--------------------------------------------------------------------------
pE  = spm_vec(DCM.M.pE);
wEp = 0;
for model = 1:N
    wEp   = wEp + miCp(:,:,model)*mEp(:,model);
end
Ep  = Cp*wEp;

% project back through U
%--------------------------------------------------------------------------
Cp  = U*Cp*U';
Ep  = U*Ep + pE - U*U'*pE;
Ep  = spm_unvec(Ep,DCM.M.pE);


%-Copy contents of first DCM into the output DCM and add BPA
%==========================================================================
DCM            = DCM_first;
DCM.averaged   = true;
try
    DCM.models = char(P);
end

% compute posterior probabilities and variance
%--------------------------------------------------------------------------
sw      = warning('off','SPM:negativeVariance');
Vp      = diag(Cp);
Pp      = 1 - spm_Ncdf(0,abs(spm_vec(Ep) - spm_vec(pE)),Vp);
warning(sw);

DCM.Ep  = Ep;
DCM.Cp  = Cp;
DCM.Vp  = spm_unvec(Vp,DCM.M.pE);
DCM.Pp  = spm_unvec(Pp,DCM.M.pE);


%-Save new DCM if there are no output arguments
%==========================================================================
if nocond
    DCM.name = [name ' (FFX average - no conditional dependencies)'];
else
    DCM.name = [name ' (Bayesian FFX average)'];
end

if ~nargout
    save(['DCM_avg_' name '.mat'], 'DCM', spm_get_defaults('mat.format'));
end

% Warn the user how this average DCM should NOT be used
%--------------------------------------------------------------------------
% disp(['Results of averaging DCMs were saved in DCM_avg_' name '.mat.']);
% disp('Please note that this file only contains average parameter estimates');
% disp('and their posterior probabilities, but NOT averaged time series.');
% disp('Also, note that this file can NOT be used for model comparisons.');

