function DCM = spm_dcm_erp_data(DCM)
% prepares structures for forward model (both EEG and MEG)
% FORMAT DCM = spm_dcm_erp_data(DCM)
% DCM   - DCM structure
% requires
%
%    DCM.xY.Dfile
%    DCM.options.trials
%    DCM.options.Tdcm
%    DCM.options.D
%
% sets
%    DCM.xY.Time - Time [ms] of downsampled data
%    DCM.xY.dt   - sampling in seconds
%    DCM.xY.y    - concatenated response
%    DCM.xY.It   - Indices of time bins
%    DCM.xY.Ic   - Indices of good channels
%__________________________________________________________________________
% Stefan Kiebel, Karl friston
% $Id: spm_dcm_erp_data.m 668 2006-10-26 16:35:28Z karl $

% Get D filename
%--------------------------------------------------------------------------
try
    Dfile = DCM.xY.Dfile;
catch
    errordlg('Please specify data and trials');
    error('')
end

% load D
%--------------------------------------------------------------------------
D = spm_eeg_ldata(Dfile);

% indices of EEG channel (excluding bad channels) and perstimulus times
%--------------------------------------------------------------------------
Ic              = setdiff(D.channels.eeg, D.channels.Bad);
DCM.M.dipfit.Ic = Ic;
DCM.xY.Ic       = Ic;
DCM.xY.Time     = 1000*[-D.events.start:D.events.stop]/D.Radc; % ms
DCM.xY.dt       = 1/D.Radc;

% options
%--------------------------------------------------------------------------
try
    DT   = DCM.options.D;
catch
    errordlg('Please specify down sampling');
    error('')
end
try
    T1   = DCM.options.Tdcm(1);
    T2   = DCM.options.Tdcm(2);
catch
    errordlg('Please specify time window');
    error('')
end


% if MEG, store grad struct in D.channels
%--------------------------------------------------------------------------
try
    DCM.xY.grad = D.channels.grad;
end

% time window and bins for modelling
%--------------------------------------------------------------------------
[i, T1] = min(abs(DCM.xY.Time - T1));
[i, T2] = min(abs(DCM.xY.Time - T2));
It      = [T1:DT:T2]';                    % time bins
try
    trial     = DCM.options.trials;
    DCM.xY.xy = {};
    for i = 1:length(trial);
        DCM.xY.xy{i} = D.data(Ic,It,trial(i))';
    end
catch
    errordlg('please specify trials');
    error('')
end

% concatenate
%--------------------------------------------------------------------------
DCM.xY.Time = DCM.xY.Time(It);           % Time [ms] of downsampled data
DCM.xY.dt   = DT/D.Radc;                 % sampling in seconds
DCM.xY.y    = spm_cat(DCM.xY.xy(:));     % concatenated response
DCM.xY.It   = It;                        % Indices of time bins
DCM.xY.Ic   = Ic;                        % Indices of good channels

