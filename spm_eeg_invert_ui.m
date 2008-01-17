function [D] = spm_eeg_invert_ui(varargin)
% GUI for ReML inversion of forward model for EEG-EMG
% FORMAT [D] = spm_eeg_invert_ui(D,val)
% ReML estimation of regularisation hyperparameters using the
% spatio-temporal hierarchy implicit in EEG data
% sets:
%
%     D.inv{i}.inverse.trials - trials (in D.events.types) to invert
%     D.inv{i}.inverse.smooth - smoothness of source priors (mm)
%     D.inv{i}.inverse.type   - 'MSP' multiple sparse priors
%                               'LOR' LORETA-like model
%                               'IID' LORETA and WMN
%     D.inv{i}.inverse.xyz    - (n x 3) locations of spherical VOIs
%     D.inv{i}.inverse.rad    - radius (mm) of VOIs
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id:$

% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});


% check whether to use conventional or DCM temporal priors
%--------------------------------------------------------------------------
if spm_input('Reconstruction','+1','b',{'Classical|DCM'},[0 1],1)
    
        % record type in D and DCM structures
        %------------------------------------------------------------------
        inverse.type = 'DCM';

        % exchange filenames
        %------------------------------------------------------------------
        DCMfile            = ['DCM_' D.fname];
        D.inv{val}.DCMfile = DCMfile;
        DCM.val            = val;
        DCM.xY.Dfile       = fullfile(D.path,D.fname);
        DCM.options.type   = 2;
        DCM.name           = DCMfile;
        
        % an call API to specify DCM
        %------------------------------------------------------------------
        spm_api_erp(DCM);
        D.inv{val}.inverse = inverse;
        return
end

% Conventional reconstruction: get conditions or trials
%==========================================================================
if length(D.events.types) > 1
    if spm_input('All conditions or trials','+1','b',{'yes|no'},[1 0],1)
        trials = D.events.types;
    else
        trials = [];
        for  i = 1:length(D.events.types)
            str = sprintf('invert %i',D.events.types(i))
            if spm_input(str,'+1','b',{'yes|no'},[1 0],1);
                trials(end + 1) = D.events.types(i);
            end
        end
    end
else
    trials = D.events.types;
end

% Inversion parameters
%--------------------------------------------------------------------------
inverse        = spm_eeg_inv_custom_ui(D);
inverse.trials = trials;

% invert
%==========================================================================
D.con               = 1;
D.inv{val}.inverse  = inverse;
D                   = spm_eeg_invert(D);
