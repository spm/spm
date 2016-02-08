function [P]   = spm_dcm_fit(P)
% Bayesian inversion of DCMs using Variational Laplace
% FORMAT [DCM] = spm_dcm_fit(P)
%
% P    - {N x M} DCM structure array (or filenames) from N subjects
%
% DCM  - Inverted (1st level) DCM structures with posterior densities
%__________________________________________________________________________
%
% This routine is just a wrapper that calls the appropriate dcm inversion
% routine for a set a pre-specifed DCMs.
%
% If called with a cell array, each column is assumed to contain 1st level
% DCMs inverted under the same model. Each row contains a different data
% set (or subject).
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fit.m 6716 2016-02-08 18:21:37Z peter $


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

model = spm_dcm_identify(DCM);

if isempty(model)    
    warning('unknown inversion scheme');
    return   
end

% get data structure for each subject (column)
%------------------------------------------------------------------
for i = 1:Ns
    for j = 2:Nm
        switch model
            
            case{'DEM'}
                P{i, j}.xY = P{i, 1}.Y;
            otherwise
                P{i, j}.xY = P{i, 1}.xY;
        end
    end
end

%matlabpool open 7

% loop over subjects (columns)
%--------------------------------------------------------------------------
%parfor i = 1:numel(P)

for i = 1:numel(P)
    
    % loop over models (rows)
    %----------------------------------------------------------------------
    
    % Get model specification
    %------------------------------------------------------------------
    try, load(P{i}); catch, DCM = P{i}; end
    
    
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
                
                % behavioural Markov decision process model
                %----------------------------------------------------------
            case{'MDP'}
                DCM = spm_dcm_mdp(DCM);
                
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
            try
                DCM = feval(model, DCM);
            catch
                error('unknown DCM');
            end
    end
    
    % place inverted model in output array
    %------------------------------------------------------------------
    P{i} = DCM;
end


%matlabpool close