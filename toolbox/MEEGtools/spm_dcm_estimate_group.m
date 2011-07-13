function spm_dcm_estimate_group(DCMs, DD, P, pE, pC)
% Apply a set of pre-specified DCMs to a set of subjects.
%
% FORMAT spm_dcm_estimate_group(DCM, D, P, pE, pC)
%  
% Arguments (optional)
%
% DCMs - a list of DCM files
% DD   - a list of MEEG datasets
% P   - initialisation (1 - use previous posteriors)
% pE  - priors (1 - take from DCM)
% pC  - prior covariance
%
% All results will be saved in the current directory
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_dcm_estimate_group.m 4390 2011-07-13 18:04:22Z vladimir $

% Disclaimer: this code is provided as an example and is not guaranteed to
% work with data on which it was not tested. If it does not work for you,
% feel free to improve it and contribute your improvements to the MEEGtools
% toolbox in SPM (http://www.fil.ion.ucl.ac.uk/spm/)

if nargin == 0
    DCMs =  spm_select(Inf, 'mat', 'Select DCM mat files');
end

if nargin < 2
    DD =  spm_select(Inf, 'mat', 'Select SPM M/EEG mat files');
end

if nargin < 3,    P  = [];  end
if nargin < 4,    pE = [];  end
if nargin < 5,    pC = [];  end

for i = 1:size(DCMs, 1)
    cDCM = getfield(load(deblank(DCMs(i, :)), 'DCM'), 'DCM');
    
    % initialise with posteriors if required
    % -------------------------------------------------------------------------
    if isequal(P, 1)
        cDCM.M.P = cDCM.Ep;
    else
        cDCM.M.P = P;
    end
    
    % initialise with posteriors if required
    % -------------------------------------------------------------------------
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
        
        [ok, D] = check(D, 'dcm');
        
        if ~ok
            if check(D, 'basic')
                warning (['The file ' D.fname ' is not ready for DCM.'...
                    'Use prepare to specify sensors and fiducials or LFP channels.']);
            else
                warning(['The meeg file ' D.fname ' is corrupt or incomplete']);
            end
            continue;
        end
        
        
        DCM.xY.Dfile = fullfile(D.path, D.fname);
        
        DCM  = spm_dcm_erp_data(DCM, DCM.options.h);
        
        [p, f] = fileparts(DCM.name);
        
        DCM.name = fullfile(pwd, [f '_' D.fname]);
        
        % invert and save
        %--------------------------------------------------------------------------
        switch DCM.options.analysis
            
            % conventional neural-mass and mean-field models
            %----------------------------------------------------------------------
            case{'ERP'}
                DCM = spm_dcm_erp(DCM);
                
                % cross-spectral density model (complex)
                %----------------------------------------------------------------------
            case{'CSD'}
                DCM = spm_dcm_csd(DCM);
                
                % cross-spectral density model (steady-state responses)
                %----------------------------------------------------------------------
            case{'SSR'}
                DCM = spm_dcm_ssr(DCM);
                
                % induced responses
                %----------------------------------------------------------------------
            case{'IND'}
                DCM = spm_dcm_ind(DCM);
                
                % phase coupling
                %----------------------------------------------------------------------
            case{'PHA'}
                DCM = spm_dcm_phase(DCM);
                
            otherwise
                error('unknown DCM type')
        end
    end
end
        