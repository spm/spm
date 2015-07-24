function [P]   = spm_dcm_fit(P)
% Bayesian inversion of DCMs using Variational Laplace
% FORMAT [DCM] = spm_dcm_fit(P)
%
% P    - {N x M} DCM structure array (or filenames) from N subjects
%
% DCM  - Inverted (1st level) DCM structures with posterior densities
%__________________________________________________________________________
%
% This routine is just a wrapper that calls the apprioriate dcm inversion
% routine for a set a pre-specifed DCMs.
%
% If called with a cell array, each column is assumed to contain 1st level
% DCMs inverted under the same model. Each row contains a different data
% set (or subject).
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fit.m 6506 2015-07-24 10:26:51Z karl $


% get filenames and set up
%--------------------------------------------------------------------------
if ~nargin
    [P, sts] = spm_select([1 Inf],'^DCM.*\.mat$','Select DCM*.mat files');
    if ~sts, return; end
end
if ischar(P),   P = cellstr(P);  end
if isstruct(P), P = {P};         end

% Number of subjects (data) and models (of those data)
%--------------------------------------------------------------------------
[Ns,Nm] = size(P);

% Find model class and modality
%==========================================================================
try, load(P{1}); catch, DCM = P{1}; end
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
    
elseif isfield(DCM.M,'IS')
    
    % assume the model is specified explicitly
    %----------------------------------------------------------------------
    model  = 'NLSI';
    
elseif isfield(DCM.M,'E')
    
    % assume the model is a hierarchical dynamic model
    %----------------------------------------------------------------------
    model  = 'DEM';
    
else
    
    warning('unknown inversion scheme');
    return
    
end

% loop over subjects (columns)
%--------------------------------------------------------------------------
for i = 1:Ns
    
    % loop over models (rows)
    %----------------------------------------------------------------------
    for j = 1:Nm
        
        % Get model specification
        %------------------------------------------------------------------
        try, load(P{i,j}); catch, DCM = P{i,j}; end
        
        % get data structure for this subject (column)
        %------------------------------------------------------------------
        switch model
            
            case{'DEM'}
                if j == 1, Y  = DCM.Y;  else, DCM.Y  = Y;  end
            otherwise
                if j == 1, xY = DCM.xY; else, DCM.xY = xY; end
                
        end
        
        % invert and save
        %==================================================================
        switch model
            
            % fMRI model
            %--------------------------------------------------------------
            case{'fMRI'}
                DCM = spm_dcm_estimate(DCM);
                
                % conventional neural-mass and mean-field models
                %----------------------------------------------------------
            case{'fMRI_CSD'}
                DCM = spm_dcm_fmri_csd(DCM);
                
                % conventional neural-mass and mean-field models
                %----------------------------------------------------------
            case{'ERP'}
                DCM = spm_dcm_erp(DCM);
                
                % cross-spectral density model (complex)
                %----------------------------------------------------------
            case{'CSD'}
                DCM = spm_dcm_csd(DCM);
                
                % cross-spectral density model (steady-state responses)
                %----------------------------------------------------------
            case{'TFM'}
                DCM = spm_dcm_tfm(DCM);
                
                % induced responses
                %----------------------------------------------------------
            case{'IND'}
                DCM = spm_dcm_ind(DCM);
                
                % phase coupling
                %----------------------------------------------------------
            case{'PHA'}
                DCM = spm_dcm_phase(DCM);
                
                % cross-spectral density model (steady-state responses)
                %----------------------------------------------------------
            case{'NFM'}
                DCM = spm_dcm_nfm(DCM);
                
                % generic nonlinear system identification
                %----------------------------------------------------------
            case{'NLSI'}
                [Ep,Cp,Eh,F] = spm_nlsi_GN(DCM.M,DCM.xU,DCM.xY);
                DCM.Ep       = Ep;
                DCM.Eh       = Eh;
                DCM.Cp       = Cp;
                DCM.F        = F;
                
                % hierarchical ddynamic mmodel
                %----------------------------------------------------------
            case{'DEM'}
                DCM = spm_DEM(DCM);

                
            otherwise
                spm('alert!','unknown DCM','Warning');
                return
        end
        
        % place inverted model in output array
        %------------------------------------------------------------------
        P{i,j} = DCM;
        
    end
end
