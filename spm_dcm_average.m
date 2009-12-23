function spm_dcm_average (mtype,P,name,Nowarning)
% Produce an aggregate DCM model using Bayesian averaging
% FORMAT spm_dcm_average (mtype,P,name,Nowarning)
%
% mtype        -  ERP: mtype =0;  fMRI: mtype > 0
% P            -  Array of DCM filenames e.g., P(1,:)='DCM1', P(2,:)='DCM2'
% name         -  Name of DCM output file. This is prefixed by 'DCM_avg_'.
% Nowarning    -  Send warning to user (default) or not
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
% A Bayesian random effects analysis can be implemented for a
% particular contrast using the spm_dcm_sessions.m function
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Will Penny & Klaas Enno Stephan
% $Id: spm_dcm_average.m 3654 2009-12-23 20:09:54Z karl $
 
 
if nargin <= 1
    % Function called without parameters (e.g. via GUI)
    %----------------------------------------------------------------------
    Finter = spm_figure('GetWin','Interactive');
    set(Finter,'name','Dynamic Causal Modeling')
    header = get(Finter,'Name');
    WS     = spm('WinScale');
    num_models = spm_input('How many DCM models to average ? ','+1','r',[],1);
    P      = spm_select(num_models,'^DCM.*\.mat$','Select DCM*.mat files');
    name   = spm_input('Name for DCM_avg_???.mat','+1','s');
    if nargin < 1
        mtype  = spm_input('fMRI or ERP?',1,'b',{'fMRI','ERP'},[1 0]);
    end
else
    num_models = size(P,1);
end
 
if nargin<4
    Nowarning = 0;
end
 
% Loop through all selected models and get posterior means and precisions
%--------------------------------------------------------------------------
for model = 1:num_models,
    load(P(model,:),'-mat');
 
    % Only look at those parameters with non-zero prior variance
    %----------------------------------------------------------------------
    pCdiag = diag(DCM.M.pC);
    wsel   = find(pCdiag);
 
    if model == 1
        wsel_first = wsel;
        DCM_first  = DCM;
    else
        if ~(length(wsel) == length(wsel_first))
            disp('Error in spm_dcm_average: DCMs must have same structure');
            return
        end
        if ~(wsel == wsel_first)
            disp('Error in spm_dcm_average: DCMs must have same structure');
            return
        end
    end
 
    % Get posterior precision matrix and mean
    %-------------------------------------------------------------------
    Cp              = DCM.Cp;
    Ep              = spm_vec(DCM.Ep);
    miCp(:,:,model) = inv(full(Cp(wsel,wsel)));
    mEp(:,model)    = full(Ep(wsel));
 
end
 
 
% Average models using Bayesian fixed-effects analysis -> average Ep,Cp
%==========================================================================
 
% averaged posterior covariance
%--------------------------------------------------------------------------
Cp(wsel,wsel) = inv(sum(miCp,3));
 
% averaged posterior mean
%--------------------------------------------------------------------------
pE          = DCM.pE;
weighted_Ep = 0;
for model = 1:num_models,
    weighted_Ep = weighted_Ep + miCp(:,:,model)*mEp(:,model);
end
Ep(wsel)    = Cp(wsel,wsel)*weighted_Ep;
Ep          = spm_unvec(Ep,pE);
 
% Copy contents of first DCM into the output DCM and add BPA
%==========================================================================
DCM          = DCM_first;
DCM.models   = P;
DCM.averaged = 1;
 
% compute posterior probabilities sand variance
%--------------------------------------------------------------------------
sw      = warning('off','SPM:negativeVariance');
Vp      = diag(Cp);
Pp      = 1 - spm_Ncdf(0,abs(spm_vec(Ep) - spm_vec(pE)),Vp);
warning(sw);
 
DCM.Ep  = Ep;
DCM.Cp  = Cp;
DCM.Vp  = spm_unvec(Vp,pE);
DCM.Pp  = spm_unvec(Pp,pE);
 
 
% Save new DCM
%--------------------------------------------------------------------------
DCM.name = [name ' (Bayesian FFX average)'];
if spm_matlab_version_chk('7') >= 0
    save(['DCM_avg_' name], 'DCM', '-V6');
else
    save(['DCM_avg_' name], 'DCM');
end;
 
% If called through GUI, prompt user that averaging is finished
%--------------------------------------------------------------------------
if nargin <= 1
    spm_clf
    spm('FigName',header);
    spm('Pointer','Arrow')
    spm_input(['Results of averaging DCMs were saved in DCM_avg_' name],1,'d')
end
 
% Warn the user how this average DCM should NOT be used
%--------------------------------------------------------------------------
if ~Nowarning
    str = {['Results of averaging DCMs were saved in DCM_avg_' name '.'], ...
            ' ', ...
            'Please note that this file only contains average parameter estimates and their post. probabilities, but NOT averaged time series.', ...
            ' ', ...
            'Also, note that this file can NOT be used for model comparisons.'};
    spm('alert!',str,'DCM average warning')
end
