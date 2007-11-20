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
inverse.trials = trials;


% Type of analysis
%--------------------------------------------------------------------------
type           = spm_input('Type of inversion','+1','MSP|COH|MNM',{'GS','LOR','IID'},1);
inverse.type   = type{:};

if spm_input('Model','+1','b',{'Standard|Custom'},[0 1],1)


    % D.inverse.smooth - smoothness of source priors (mm)
    %----------------------------------------------------------------------
    switch inverse.type, case{'GS','MSP','LOR'}
        inverse.smooth = spm_input('Spatial smoothness (0-1)','+1','0.2|0.4|0.6',[0.2 0.4 0.6],3);
    end
    
    % D.inverse.sdv    - smoothness of source priors (ms)
    %----------------------------------------------------------------------
    inverse.sdv      = spm_input('Temporal smoothness (ms)','+1','1|4|16',[1 4 16],2);

    % Search strategy
    %----------------------------------------------------------------------
    switch inverse.type, case{'MSP','GS'}
        type         = spm_input('Search','+1','Sparse|Greedy',{'MSP','GS'},2);
        inverse.type = type{:};
    end
    
    % Number of sparse priors
    %----------------------------------------------------------------------
    switch inverse.type, case{'MSP','GS'}
        inverse.Np   = spm_input('MSPs per hemisphere','+1','64|128|256|512',[64 128 256 512],3);
    end
    
    % Time window of interest
    %----------------------------------------------------------------------
    woi         = round([-D.events.start D.events.stop]*1000/D.Radc);
    woi         = spm_input('Time window (ms)','+1','r',woi);
    inverse.woi = round([min(woi) max(woi)]);
    
    % High-pass filter
    %----------------------------------------------------------------------
    inverse.lpf = spm_input('High-pass (Hz)','+1','1|8|16',[1 8 16],1);
    
    % Low-pass filter
    %----------------------------------------------------------------------
    inverse.hpf = spm_input('Low-pass (Hz)','+1','64|128|256',[64 128 256],3);
        
    % Source space restictions
    %----------------------------------------------------------------------
    if spm_input('Restrict solutions','+1','yes|no',[1 0],2);

        [f,p]       = uigetfile('*.mat','source (n x 3) location file');
        xyz         = load(fullfile(p,f));
        name        = fieldnames(xyz);
        xyz         = getfield(xyz, name{1});
        inverse.xyz = xyz;
        inverse.rad = spm_input('radius of VOI (mm)','+1','r',32);
    end

end

% invert
%==========================================================================
D.con               = 1;
D.inv{val}.inverse  = inverse;
D                   = spm_eeg_invert(D);


return
%==========================================================================
% other GUI options

    % High-pass filter
    %----------------------------------------------------------------------
    inverse.Han = spm_input('PST Hanning','+1','yes|no',[1 0],1);

    % Low-pass filter
    %----------------------------------------------------------------------
    inverse.hpf = spm_input('Low-pass (Hz)','+1','64|128|256',[64 128 256],3);
        
    % High-pass filter
    %----------------------------------------------------------------------
    inverse.lpf = spm_input('High-pass (Hz)','+1','1|8|32',[1 8 32],1);

    % Channel modes
    %----------------------------------------------------------------------
    inverse.Nm   = spm_input('Channel modes (max)','+1','32|64|128',[32 64 128],2);


