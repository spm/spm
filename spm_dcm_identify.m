function model = spm_dcm_identify(DCM)
% Identify the type of DCM. Return an empty string if unknown.
%
% DCM   - the model to evaluate
%
% model - a string identifying the modality
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging
%
% $Id: spm_dcm_identify.m 6716 2016-02-08 18:21:37Z peter $

if ischar(DCM)
    DCM = load(DCM);
    DCM = DCM.DCM;
end

if isfield(DCM,'options')
    
    % spatial models for EEG
    %----------------------------------------------------------------------
    if isfield(DCM.options,'spatial')
        model = DCM.options.analysis;
    else
        
        % an fMRI model
        %------------------------------------------------------------------
        try
            DCM.options.analysis;
            model = 'fMRI_CSD';
        catch
            model = 'fMRI';
        end
        
    end
    
elseif isfield(DCM,'MDP')
    
    % assume the model is specified explicitly
    %----------------------------------------------------------------------
    model  = 'MDP';
    
elseif isfield(DCM.M,'IS')
    
    % assume the model is specified explicitly
    %----------------------------------------------------------------------
    model  = 'NLSI';
    
elseif isfield(DCM.M,'E')
    
    % assume the model is a hierarchical dynamic model
    %----------------------------------------------------------------------
    model  = 'DEM';
    
else
    
    model = '';
    return
    
end