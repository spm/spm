function [] = spm_dcm_average (num_models,P,name)
% Produce an aggregate DCM model using Bayesian averaging
% FORMAT [] = spm_dcm_average (num_models,P,name)
%
% num_models        The number of models
% P                 Cell array of DCM filenames eg. P{1}='DCM1', P{2}='DCM2'
% name              Name of DCM output file. This is prefixed by 'DCM_avg_'.
%
% This routine creates a new DCM model which is the average over a 
% number of fitted DCM models. These can be over sessions or over subjects.
% This average model can then be interrogated using the standard 
% DCM 'review' options eg. to look at individual parameters or 
% contrasts of parameters. The resulting inferences correspond to 
% a Bayesian Fixed Effects analysis.
%
% Note that the Bayesian averaging is only applied to the A, B and C matrices.
% All other quantities in the average model are simply copied from 
% the first DCM in the list. Only models with exactly the same 
% A,B,C structure can be averaged.
%
% A Bayesian random effects analysis can be implemented for a 
% particular contrast using the spm_dcm_sessions.m function
%
% -------------------------------------------------------
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Will Penny
% $Id$


if nargin < 1
    % Function called via UI
    
    %-display model details
    Finter = spm_figure('GetWin','Interactive');
    %-------------------------------------------------------------------
    set(Finter,'name','Dynamic Causal Modeling')
    %-get results
    %-------------------------------------------------------------------
    num_models = spm_input('How many DCM models to average ? ','+1','r',[],1);
    P     = spm_get(num_models,'DCM*.mat',{'select DCM*.mat files'});
    name  = spm_input('name for DCM_avg_???.mat','+1','s');
end

% Average models using Bayesian fixed effects analysis by computing new Ep,Cp
for model=1:num_models,
    load(P{model});
    
    % Only look at those parameters with non-zero prior covariance
    pCdiag=diag(DCM.M.pC);
    wsel=find(~(pCdiag==0));
    
    if model==1
        wsel_first=wsel;
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
    DCM_first=DCM;
    
    % Only look at A,B,C values - ignore hemodynamics
    m=DCM.M.m; % Number of inputs
    n=DCM.n; % Number of regions
    % Ignore last 5*n parameters (hemodynamic coeffs)
    cwsel=wsel(1:end-5*n);
    
    % Get posterior precision matrix from model
    miCp(:,:,model)=inv(full(DCM.Cp(cwsel,cwsel)));
    % Get posterior mean from model
    mEp(:,model)=full(DCM.Ep(cwsel));
end

% Now set up average DCM model
DCM=DCM_first;
final_iCp=sum(miCp,3);
Cp=inv(final_iCp);

weighted_Ep=zeros(length(cwsel),1);
for model=1:num_models,
    weighted_Ep=weighted_Ep+miCp(:,:,model)*mEp(:,model);
end
Ep=Cp*weighted_Ep;

DCM.Ep=DCM_first.Ep;
DCM.Cp=DCM_first.Cp;
DCM.Ep(cwsel)=Ep;
DCM.Cp(cwsel,cwsel)=Cp;

% Now reshape into parameters, variances and probabilities
T          = log(2)/4;			
pp         = 1 - spm_Ncdf(T,abs(DCM.Ep),diag(DCM.Cp));
[ A  B  C] = spm_dcm_reshape(DCM.Ep,m,n,1);
[pA pB pC] = spm_dcm_reshape(pp,m,n,1);
vv         = diag(DCM.Cp);
[vA vB vC] = spm_dcm_reshape(vv,m,n,1);

% Put in data structure
DCM.A=A;
DCM.B=B;
DCM.C=C;
DCM.pA=pA;
DCM.pB=pB;
DCM.pC=pC;
DCM.vA=vA;
DCM.vB=vB;
DCM.vC=vC;

% Save new DCM
if str2num(version('-release'))>=14,
    save(['DCM_avg_' name], 'DCM', '-V6');
else
    save(['DCM_avg_' name], 'DCM');
end;

