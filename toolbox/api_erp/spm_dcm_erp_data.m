function DCM = spm_dcm_erp_data(DCM,h)
% prepares structures for forward model (both EEG and MEG)
% FORMAT DCM = spm_dcm_erp_data(DCM)
% DCM    -  DCM structure
% ERP    - 'ERP' or 'Induced' for evoked or induced response 
% requires
%
%    DCM.xY.Dfile        - data file
%    DCM.options.trials  - trial codes
%    DCM.options.Tdcm    - Peri-stimulus time window
%    DCM.options.D       - Down-sampling
%    DCM.options.han     - hanning
%    
% sets
%    DCM.xY.modality - 'MEG','EEG' or 'LFP'
%    DCM.xY.Time     - Time [ms] data
%    DCM.xY.pst      - Time [ms] of down-sampled data
%    DCM.xY.dt       - sampling in seconds (s)
%    DCM.xY.y        - response variable for DCM
%    DCM.xY.xy       - cell array of trial-speficic response {[ns x nc]}
%    DCM.xY.It       - Indices of (ns) time bins
%    DCM.xY.Ic       - Indices of (nc) good channels
%    DCM.xY.name     - names of(nc) channels

%
%    DCM.xY.Hz       - Frequency bins (for Wavelet transform)
%    DCM.options.h
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_erp_data.m 1208 2008-03-13 20:59:12Z karl $
 
% Set defaults and Get D filename
%--------------------------------------------------------------------------
try
    Dfile = DCM.xY.Dfile;
catch
    errordlg('Please specify data and trials');
    error('')
end
 
% order of drift terms
%--------------------------------------------------------------------------
try h; catch h = 0; end
 
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
        try
            [f,p]        = uigetfile('*.mat','please select data file'); 
            name         = fullfile(p,f);
            D            = spm_eeg_ldata(name);
            DCM.xY.Dfile = fullfile(name);
        catch
            warndlg([Dfile ' could not be found'])
        return
        end
    end
end
 
 
% indices of EEG channel (excluding bad channels) and peristimulus times
%--------------------------------------------------------------------------
Ic              = setdiff(D.channels.eeg, D.channels.Bad);
Nc              = length(Ic);
DCM.xY.name     = D.channels.name(Ic);
DCM.xY.modality = D.modality;
DCM.xY.Ic       = Ic;
DCM.xY.Time     = 1000*[-D.events.start:D.events.stop]/D.Radc; % PST (ms)
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
    han  = DCM.options.han;
catch
    han  = 0;
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
    T1          = DCM.options.Tdcm(1);
    T2          = DCM.options.Tdcm(2);
    [i, T1]     = min(abs(DCM.xY.Time - T1));
    [i, T2]     = min(abs(DCM.xY.Time - T2));
    
    % Time [ms] of down-sampled data
    %----------------------------------------------------------------------
    It          = [T1:DT:T2]';
    Ns          = length(It);                % number of samples
    DCM.xY.pst  = DCM.xY.Time(It);           % Down-sampled PST
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
        c = find(D.events.code == D.events.types(trial(i)) & ~D.events.reject);
    else
        c = find(D.events.code == D.events.types(trial(i)));
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


% confounds - DCT:
%--------------------------------------------------------------------------
if h == 0
    X0 = sparse(Ns,1);
else
    X0 = spm_dctmtx(Ns,h);
end
R      = speye(Ns) - X0*X0';
 
% hanning
%--------------------------------------------------------------------------
if han
    R  = R*diag(hanning(Ns))*R;
end
 
% adjust data
%--------------------------------------------------------------------------
for i = 1:length(DCM.xY.xy);
  DCM.xY.xy{i} = R*DCM.xY.xy{i};
end
 
% condition units of measurement
%--------------------------------------------------------------------------
DCM.xY.y    = spm_cond_units(DCM.xY.xy);
DCM.xY.code = D.events.code(trial);
