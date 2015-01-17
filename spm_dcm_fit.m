function [P]   = spm_dcm_fit(P)
% Bayesian inversion of DCMs using Variational Laplace
% FORMAT [DCM] = spm_dcm_fit(P)
%
% P    - {N x M} DCM structure array (or filenames) from N subjects
%
% DCM  - Inverted (1st level) DCM structures with emprical priors
%
%--------------------------------------------------------------------------
% This routine is just a wrapper that calls the apprioriate dcm inversion
% routine. for a set a pre-specifed DCMs.
%
% If called with a cell array, each column is assumed to contain 1st level
% DCMs inverted under the same model. Each row contains a different data
% set (or subject).
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fit.m 6305 2015-01-17 12:40:51Z karl $


% get filenames and set up
%--------------------------------------------------------------------------
if ~nargin
    [P, sts] = spm_select([2 Inf],'^DCM.*\.mat$','Select DCM*.mat files');
    if ~sts, return; end
end
if ischar(P),   P = cellstr(P);  end
if isstruct(P), P = {P};         end

% Number of subjects (data) and models (of those data)
%--------------------------------------------------------------------------
[Ns,Nm] = size(P);

% Find model class and modality
%--------------------------------------------------------------------------
try, load(P{1}); catch, DCM = P{1}; end
if isfield(DCM.options,'spatial')
    model = DCM.options.analysis;
else
    try
        DCM.options.analysis;
        model = 'fMRI_CSD';
    catch
        model = 'fMRI';
    end
    
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
        if j == 1
            xY = DCM.xY;
        else
            DCM.xY = xY;
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
                
            otherwise
                warndlg('unknown DCM')
                return
        end
        
        % place inverted model in output array
        %------------------------------------------------------------------
        P{i,j} = DCM;
        
    end
end
