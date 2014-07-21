function spm_dcm_estimate_group(DCMs, DD, P, pE, pC, feedback)
% Apply a set of pre-specified DCMs to a set of subjects
%
% FORMAT spm_dcm_estimate_group(DCM, D, P, pE, pC, feedback)
%  
% Arguments (optional)
%
% DCMs     - a list of DCM files
% DD       - a list of MEEG datasets
% P        - initialisation (1 - use previous posteriors)
% pE       - priors (1 - take from DCM)
% pC       - prior covariance
% feedback - provide graphical feedback [Default: 1]
%
% All results will be saved in the current directory.
%__________________________________________________________________________
% Copyright (C) 2011-2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_dcm_estimate_group.m 6112 2014-07-21 09:39:53Z karl $


if nargin == 0
    DCMs = spm_select(Inf, 'mat', 'Select DCM mat files');
end

if nargin < 2
    DD = spm_select(Inf, 'mat', 'Select SPM M/EEG mat files');
end

if nargin < 3, P  = []; end
if nargin < 4, pE = []; end
if nargin < 5, pC = []; end

if nargin < 6, feedback = 1; end


for i = 1:size(DCMs, 1)
    cDCM = getfield(load(deblank(DCMs(i, :)), 'DCM'), 'DCM');
    
    % initialise with posteriors if required
    % ---------------------------------------------------------------------
    if isequal(P, 1)
        cDCM.M.P = cDCM.Ep;
    else
        cDCM.M.P = P;
    end
    
    % initialise with posteriors if required
    % ---------------------------------------------------------------------
    if isempty(pE)
        if isfield(cDCM.M,'pE')
            cDCM.M = rmfield(cDCM.M,'pE');
        end
        if isfield(cDCM.M,'pC')
            cDCM.M = rmfield(cDCM.M,'pC');
        end
    elseif ~isequal(pE, 1)
        cDCM.M.pE = pE;
        if ~isempty(pC)
            cDCM.M.pC = pC;
        end
    end
    
    for j = 1:size(DD, 1)
        DCM = cDCM;
        
        D = spm_eeg_load(deblank(DD(j, :)));
        
        [D, ok] = check(D, 'dcm');
        
        if ~ok
            if check(D, 'basic')
                warning (['The file ' D.fname ' is not ready for DCM.'...
                    'Use prepare to specify sensors and fiducials or LFP channels.']);
            else
                warning(['The meeg file ' D.fname ' is corrupt or incomplete']);
            end
            continue;
        end
        
        
        DCM.xY.Dfile  = fullfile(D.path, D.fname);
        DCM           = spm_dcm_erp_data(DCM);
        [p, f]        = fileparts(DCM.name);
        DCM.name      = fullfile(pwd, [f '_' D.fname]);
        DCM.M.nograph = ~feedback;
        
        % invert and save
        %------------------------------------------------------------------
        switch DCM.options.analysis
            
            % conventional neural-mass and mean-field models
            %--------------------------------------------------------------
            case{'ERP'}
                DCM = spm_dcm_erp(DCM);
                
            % cross-spectral density model (complex)
            %--------------------------------------------------------------
            case{'CSD'}
                DCM = spm_dcm_csd(DCM);
                
            % cross-spectral density model (steady-state responses)
            %--------------------------------------------------------------
            case{'SSR'}
                DCM = spm_dcm_ssr(DCM);
                
            % induced responses
            %--------------------------------------------------------------
            case{'IND'}
                DCM = spm_dcm_ind(DCM);
                
            % phase coupling
            %--------------------------------------------------------------
            case{'PHA'}
                DCM = spm_dcm_phase(DCM);
                
            otherwise
                error('unknown DCM type');
        end
    end
end
        