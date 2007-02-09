function [D] = spm_eeg_invert_ui(varargin)
% GUI for ReML inversion of forward model for EEG-EMG
% FORMAT [D] = spm_eeg_invert_ui(D,val)
% ReML estimation of regularisation hyperparameters using the
% spatio-temporal hierarchy implicit in EEG data
% sets:
%
%     D.inv{i}.inverse.trials - indces of D.events.types to invert
%     D.inv{i}.inverse.con    - condition or trial type
%     D.inv{i}.inverse.smooth - smoothness of source priors (mm)
%     D.inv{i}.inverse.type   - 'MSP' multiple sparse priors
%                               'LOR' LORETA-like model
%                               'IID' LORETA and WMN
%     D.inv{i}.inverse.xyz    - (n x 3) locations of spherical VOIs
%     D.inv{i}.inverse.rad    - radius (mm) of VOIs
%__________________________________________________________________________


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
        DCM.options.type   = 3;
        DCM.name           = DCMfile;
        
        % an call API to specify DCM
        %------------------------------------------------------------------
        spm_api_erp(DCM);
        D.inv{val}.inverse = inverse;
        return
end

% Conventional reconstruction: get conditions or trials
%==========================================================================
if D.events.Ntypes > 1
    if spm_input('All conditions or trial','+1','b',{'yes|no'},[1 0],1)
        trials = D.events.types;
    else
        str    = sprintf('which condition[s] 1:(%i)',length(D.events.types))
        trials = spm_input(str,'+1','r');
    end
else
    trials = D.events.types(1);
end
inverse.trials = trials;
D.con          = 1;

% Type of analysis
%--------------------------------------------------------------------------
inverse.type   = spm_input('Type of inversion','+1','MSP|LOR|IID');

    
% D.inverse.smooth - smoothness of source priors (mm)
%--------------------------------------------------------------------------
switch inverse.type, case{'MSP','LOR'}
    inverse.smooth = spm_input('Smoothness (0-1)','+1','0.2|0.4|0.6',[0.2 0.4 0.6],2);
end

% Number of sparse priors
%--------------------------------------------------------------------------
switch inverse.type, case{'MSP'}
    inverse.Np  = spm_input('MSPs per hemisphere','+1','64|96|128',[64 96 128]);
end

% Source space restictions
%--------------------------------------------------------------------------
if spm_input('Restrict solutions','+1','yes|no',[1 0],2);
    
    [f,p]       = uigetfile('*.mat','source (n x 3) location file');
    xyz         = load(fullfile(p,f));
    name        = fieldnames(xyz);
    xyz         = getfield(xyz, name{1});
    inverse.xyz = xyz;
    inverse.rad = spm_input('radius of VOI (mm)','+1','r',32);
end

% invert
%==========================================================================
D.inv{val}.inverse = inverse;
D                  = spm_eeg_invert(D);





