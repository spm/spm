function DCM = spm_dcm_erp_data(DCM)
% prepares structures for forward model (both EEG and MEG)
% FORMAT DCM = spm_dcm_erp_data(DCM)
% DCM    -  DCM structure
% ERP    - 'ERP' or 'Induced' for evoked or induced response 
% requires
%
%    DCM.xY.Dfile
%    DCM.options.trials
%    DCM.options.Tdcm
%    DCM.options.D
%    
%
% sets
%    DCM.xY.modality - 'MEG' or 'EEG'
%    DCM.xY.Time     - Time [ms] of downsampled data
%    DCM.xY.dt       - sampling in seconds
%    DCM.xY.y        - concatenated response
%    DCM.xY.It       - Indices of time bins
%    DCM.xY.Ic       - Indices of good channels
%
%    DCM.xY.Hz       - Frequency bins (for Wavelet transform)
%    DCM.options.h
%__________________________________________________________________________
% Stefan Kiebel, Karl friston
% $Id: spm_dcm_erp_data.m 668 2006-10-26 16:35:28Z karl $

% Set defaults and Get D filename
%--------------------------------------------------------------------------
try
    Dfile = DCM.xY.Dfile;
catch
    errordlg('Please specify data and trials');
    error('')
end

% load D
%--------------------------------------------------------------------------
try
    D = spm_eeg_ldata(Dfile);
catch
    try
        [p,f]        = fileparts(Dfile);
        D            = spm_eeg_ldata(f);
        DCM.xY.Dfile = fullfile(pwd,f);
    catch
        warndlg([Dfile ' could not be found'])
        return
    end
end

% indices of EEG channel (excluding bad channels) and perstimulus times
%--------------------------------------------------------------------------
Ic              = setdiff(D.channels.eeg, D.channels.Bad);
Nc              = length(Ic);
DCM.xY.modality = D.modality;
DCM.xY.Ic       = Ic;
DCM.xY.Time     = 1000*[-D.events.start:D.events.stop]/D.Radc; % ms
DCM.xY.dt       = 1/D.Radc;
DCM.xY.xy       = {};

% options
%--------------------------------------------------------------------------
try
    DT   = DCM.options.D;
catch
    errordlg('Please specify down sampling');
    error('')
end
try
    trial = DCM.options.trials;
catch
    errordlg('please specify trials');
    error('')
end
try
    
    % time window and bins for modelling
    %----------------------------------------------------------------------
    T1      = DCM.options.Tdcm(1);
    T2      = DCM.options.Tdcm(2);
    [i, T1] = min(abs(DCM.xY.Time - T1));
    [i, T2] = min(abs(DCM.xY.Time - T2));
    
    % Time [ms] of downsampled data
    %----------------------------------------------------------------------
    It          = [T1:DT:T2]';
    Ns          = length(It);                % number of samples
    DCM.xY.Time = DCM.xY.Time(It);           % Down-sampled pst
    DCM.xY.dt   = DT/D.Radc;                 % sampling in seconds
    DCM.xY.It   = It;                        % Indices of time bins

catch
    errordlg('Please specify time window');
    error('')
end

% get trial averages - ERP
%--------------------------------------------------------------------------

for i = 1:length(trial);
    
    % trial indices
    %----------------------------------------------------------------------
    if isfield(D.events,'reject')
        c = find(D.events.code == D.events.types(i) & ~D.events.reject);
    else
        c = find(D.events.code == D.events.types(i));
    end
    Nt    = length(c);

    % ERP
    %----------------------------------------------------------------------
    Y     = zeros(Ns,Nc);
    for j = 1:Nt
        Y = Y + squeeze(D.data(Ic,It,c(j)))';
    end
    DCM.xY.xy{i} = Y/Nt;
end

% concatenate response and return (unless induced responses are required)
%--------------------------------------------------------------------------
DCM.xY.y    = spm_cat(DCM.xY.xy(:));
