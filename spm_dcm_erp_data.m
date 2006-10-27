function DCM = spm_dcm_erp_data(DCM,Dfile)
% prepares structures for ECD forward model (both EEG and MEG)
% FORMAT DCM = spm_dcm_erp_data(DCM,Dfile)
% Stefan Kiebel
% $Id: spm_dcm_erp_data.m 668 2006-10-26 16:35:28Z karl $

% Put data file in DCM.xY
%--------------------------------------------------------------------------
try
    DCM.xY.Dfile = Dfile;
catch
    try
        Dfile = DCM.xY.Dfile;
    catch
        try
            Dfile = DCM.M.Dfile;
            DCM.xY.Dfile = Dfile;
        catch
        [f p] = uigetfile({'*.mat'},'Please data file');
        DCM.xY.Dfile = fullfile(p,f);
        Dfile = DCM.xY.Dfile;
        end
    end
end

try
    D = spm_eeg_ldata(Dfile);
catch
    warndlg('unable to read data file')
    return
end

% indices of EEG channel (excluding bad channels) and perstimulus times
%--------------------------------------------------------------------------
chansel        = setdiff(D.channels.eeg, D.channels.Bad);
DCM.M.dipfit.chansel = chansel;
DCM.xY.chansel = chansel;
DCM.xY.Time    = 1000*[-D.events.start:D.events.stop]/D.Radc; % ms
DCM.xY.dt      = 1/D.Radc;

% options
%--------------------------------------------------------------------------
try, DT = DCM.options.D;       catch, DT  = 1;                 end
try, T1 = DCM.options.Tdcm(1); catch, T1 = DCM.xY.Time(1);    end
try, T2 = DCM.options.Tdcm(2); catch, T2 = DCM.xY.Time(end);  end


% if MEG, store grad struct in D.channels
%--------------------------------------------------------------------------
try
    DCM.xY.grad = D.channels.grad;
end

% time window and bins for modelling
%--------------------------------------------------------------------------
[i, T1] = min(abs(DCM.xY.Time - T1));
[i, T2] = min(abs(DCM.xY.Time - T2));
j       = [T1:DT:T2]';                    % time bins
try
    trial     = DCM.options.trials;
    DCM.xY.xy = {};
    for i = 1:length(trial);
        DCM.xY.xy{i} = D.data(chansel,j,trial(i))';
    end
catch
    warndlg('please specify appropriate trials');
    return
end

% concatenate
%--------------------------------------------------------------------------
DCM.xY.Time = DCM.xY.Time(j);            % Time [ms] of downsampled data
DCM.xY.dt   = DT/D.Radc;                 % sampling in seconds
DCM.xY.y    = spm_cat(DCM.xY.xy(:));     % concatenated response


