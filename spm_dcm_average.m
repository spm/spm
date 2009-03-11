function spm_dcm_average (mtype,P,name)
% Produce an aggregate DCM model using Bayesian averaging
% FORMAT spm_dcm_average (mtype,P,name)
%
% mtype        -  ERP: mtype =0;  fMRI: mtype > 0
% P            -  Array of DCM filenames eg. P(1,:)='DCM1', P(2,:)='DCM2'
% name         -  Name of DCM output file. This is prefixed by 'DCM_avg_'.
%
% This routine creates a new DCM model in which the parameters are averaged
% over a number of fitted DCM models. These can be over sessions or over 
% subjects. This average model can then be interrogated using the standard 
% DCM 'review' options to look at contrasts of parameters. The resulting 
% inferences correspond to a Bayesian Fixed Effects analysis.
%
% Note that the Bayesian averaging is only applied to the A, B and C matrices.
% All other quantities in the average model are initially simply copied from 
% the first DCM in the list. Subsequently, they are deleted before saving 
% the average DCM in order to avoid any false impression that averaged 
% models could be used for model comparison or contained averaged timeseries.
% Neither operation is valid and will be prevented by the DCM interface.
% Finally, note that only models with exactly the same A,B,C structure 
% and the same brain regions can be averaged.
%
% A Bayesian random effects analysis can be implemented for a 
% particular contrast using the spm_dcm_sessions.m function
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny & Klaas Enno Stephan
% $Id: spm_dcm_average.m 2863 2009-03-11 20:25:33Z guillaume $


if nargin <= 1
    % Function called without parameters (e.g. via GUI)
    %----------------------------------------------------------------------
    Finter = spm_figure('GetWin','Interactive');
    set(Finter,'name','Dynamic Causal Modeling')
    header = get(Finter,'Name');
    WS     = spm('WinScale');
    num_models = spm_input('How many DCM models to average ? ','+1','r',[],1);
    P     = spm_select(num_models,'^DCM.*\.mat$','Select DCM*.mat files');
    name  = spm_input('Name for DCM_avg_???.mat','+1','s');
    if nargin < 1
        mtype  = spm_input('fMRI or ERP?',1,'b',{'fMRI','ERP'},[1 0]);
    end
else
    num_models = size(P,1);
end


% Loop through all selected models and get posterior means and variances
% -------------------------------------------------------------------------
for model = 1:num_models,
    load(P(model,:),'-mat');
    
    % Only look at those parameters with non-zero prior covariance
    pCdiag = diag(DCM.M.pC);
    wsel = find(~(pCdiag==0));
    
    if model == 1
        wsel_first = wsel;
        DCM_first = DCM;
        % determine number of inputs and regions
        if mtype
            % DCM for fMRI
            m = DCM.M.m; % number of inputs
            n = DCM.n; % number of regions
            % Only look at A,B,C values - ignore hemodynamics (last 6*n parameters)
            cwsel = wsel(1:end-6*n);
        else
            % DCM for ERP
            m = size(DCM.Qp.C,2); % number of inputs
            n = size(DCM.Qp.C,1); % number of sources
            cwsel = wsel(1:end);
        end
    else
        if ~(length(wsel)==length(wsel_first))
            disp('Error in spm_dcm_average: DCM models must have same input, intrinsic and modulatory structure');
            return
        end
        if ~(wsel==wsel_first)
            disp('Error in spm_dcm_average: DCM models must have same input, intrinsic and modulatory structure');
            return
        end
    end
    
    % Get posterior precision matrix from model
    miCp(:,:,model) = inv(full(DCM.Cp(cwsel,cwsel)));
    % Get posterior mean from model
    mEp(:,model) = full(DCM.Ep(cwsel));
end


% Average models using Bayesian fixed effects analysis -> average Ep,Cp
%--------------------------------------------------------------------------
% averaged posterior covariance
final_iCp = sum(miCp,3);
Cp        = inv(final_iCp);
% averaged posterior mean
weighted_Ep = zeros(length(cwsel),1);
for model = 1:num_models,
    weighted_Ep = weighted_Ep + miCp(:,:,model)*mEp(:,model);
end
Ep = Cp*weighted_Ep;

% Copy contents of first DCM into the output DCM and insert averaged values into parameter & covariance vectors
DCM                 = DCM_first;
DCM.models          = P;
DCM.Ep(cwsel)       = Ep;
DCM.Cp(cwsel,cwsel) = Cp;

% Now reshape into parameters, variances and probabilities
if mtype
    % DCM for fMRI
    [ A  B  C] = spm_dcm_reshape(DCM.Ep,m,n,1);
    T          = 0;         
    sw = warning('off','SPM:negativeVariance');
    pp         = 1 - spm_Ncdf(T,abs(DCM.Ep),diag(DCM.Cp));
    warning(sw);
    [pA pB pC] = spm_dcm_reshape(pp,m,n,1);
    vv         = diag(DCM.Cp);
    [vA vB vC] = spm_dcm_reshape(vv,m,n,1);
    % store in DCM data structure
    DCM.A = A;
    DCM.B = B;
    DCM.C = C;
    DCM.pA = pA;
    DCM.pB = pB;
    DCM.pC = pC;
    DCM.vA = vA;
    DCM.vB = vB;
    DCM.vC = vC;    
else
    % DCM for ERP
    % get priors (note: these are the priors of the first model)
    if isfield(M, 'dipfit')
        % model with parameterised leadfield
        [pE,pC] = spm_erp_priors(DCM.A,DCM.B,DCM.C,M.dipfit,DCM.xU.dur);
    else
        % model w/ static leadfield
        L       = xY.S'*DCM.L;
        [pE,pC] = spm_erp_priors(DCM.A,DCM.B,DCM.C,M.dipfit,DCM.xU.dur);
    end
    % store in DCM data structure
    sw = warning('off','SPM:negativeVariance');
    Pp = 1 - spm_Ncdf(0,abs(spm_vec(DCM.Ep) - spm_vec(pE)),diag(DCM.Cp));
    warning(sw);
    DCM.Pp = spm_unvec(Pp,pE);
    DCM.Qp = DCM.Ep;
    DCM.averaged = 1;
end


% Remove or rename fields with subject-specific information (e.g. predicted/observed data)
% to prevent the impression that these have also been averaged and set a flag 
% indicating that this is an averaged DCM.
if mtype
    % DCM for fMRI
    try
        % DCMs computed with SPM5
        DCM.Y = rmfield (DCM.Y,{'y','X0','Ce'});
        DCM.U = rmfield (DCM.U,{'u'});
        DCM = rmfield (DCM,{'y','xY','F','AIC','BIC','R','H1','H2','K1','K2','Ce','T'});
    catch
        % minimum: remove time series
        DCM = rmfield (DCM,{'y','Y'});
    end
else
    % DCM for ERPs
    try
        DCM = rmfield (DCM,{'options','Y','U','L','xY','xU','H','Hc','K','R','Rc','Ce','F'});
    catch
        DCM = rmfield (DCM,{'Y','xY'});
    end
end
DCM.averaged = 1;
    
    
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
str = {['Results of averaging DCMs were saved in DCM_avg_' name '.'], ...
        ' ', ...
        'Please note that this file only contains average parameter estimates and their post. probabilities, but NOT averaged time series.', ...
        ' ', ...
        'Also, note that this file can NOT be used for model comparisons.'};
spm_input(str,1,'bd','OK',[1],1);
