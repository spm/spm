function BSD = spm_bsd_data(BSD)
% Cross-spectral density data-features using a VAR model
% FORMAT BSD = spm_BSD_csd_data(BSD)
% BSD    -  BSD structure
% requires
%
%    BSD.xY.Dfile        - name of data file
%    BSD.M.U             - channel subspace
%    BSD.options.trials  - trial to evaluate
%    BSD.options.Tdcm    - time limits
%    BSD.options.Fdcm    - frequency limits
%    BSD.options.D       - Down-sampling
%
% sets
%
%    BSD.xY.pst     - Peristimulus Time [ms] sampled
%    BSD.xY.dt      - sampling in seconds [s] (down-sampled)
%    BSD.xY.U       - channel subspace
%    BSD.xY.y       - cross spectral density over sources
%    BSD.xY.csd     - cross spectral density over sources
%    BSD.xY.It      - Indices of time bins
%    BSD.xY.Ic      - Indices of good channels
%    BSD.xY.Hz      - Frequency bins
%    BSD.xY.code    - trial codes evaluated
%__________________________________________________________________________
%
% Original code: spm_dcm_csd.m - Karl Friston
% Adapted by: Johan Medrano
% Copyright (C) 2024-2025 Wellcome Centre for Human Neuroimaging

 
% Set defaults and Get D filename
%-------------------------------------------------------------------------
try
    Dfile = BSD.xY.Dfile;
catch
    errordlg('Please specify data and trials');
    error('')
end

% ensure spatial modes have been computed (see spm_BSD_csd)
%-------------------------------------------------------------------------
try
    BSD.M.U;
catch
    Nm      = BSD.options.Nmodes;
    BSD.M.U = spm_dcm_eeg_channelmodes(BSD.M.dipfit,Nm);
end

% load D
%--------------------------------------------------------------------------
try
    D = spm_eeg_load(Dfile);
catch
    try
        [p,f]        = fileparts(Dfile);
        D            = spm_eeg_load(f);
        BSD.xY.Dfile = fullfile(pwd,f);
    catch
        try
            [f,p]        = uigetfile('*.mat','please select data file');
            name         = fullfile(p,f);
            D            = spm_eeg_load(name);
            BSD.xY.Dfile = fullfile(name);
        catch
            warndlg([Dfile ' could not be found'])
            return
        end
    end
end

 
% indices of EEG channel (excluding bad channels)
%--------------------------------------------------------------------------
if ~isfield(BSD.xY, 'modality')
    [mod, list] = modality(D, 0, 1);

    if isequal(mod, 'Multimodal')
        qstr = 'Only one modality can be modelled at a time. Please select.';
        if numel(list) < 4
            
            % Nice looking dialog
            %--------------------------------------------------------------
            options = [];
            options.Default = list{1};
            options.Interpreter = 'none';
            BSD.xY.modality = questdlg(qstr, 'Select modality', list{:}, options);
        else
            
            % accomodate more buttons
            %--------------------------------------------------------------
            ind = menu(qstr, list);
            BSD.xY.modality = list{ind};
        end
    else
        BSD.xY.modality = mod;
    end
end



if ~isfield(BSD.xY, 'Ic')
    BSD.xY.Ic  = D.indchantype(BSD.xY.modality,'GOOD');
end

Ic        = BSD.xY.Ic;
Nm        = size(BSD.M.U,2);
BSD.xY.Ic = Ic;

% options
%--------------------------------------------------------------------------
try
    DT    = BSD.options.D;
catch
    DT    = 1;
end
try
    trial = BSD.options.trials;
catch
    trial = 1:D.nconditions;
end
 
% check data are not oversampled (< 4ms)
%--------------------------------------------------------------------------
if DT/D.fsample < 0.004
    DT            = ceil(0.004*D.fsample);
    BSD.options.D = DT;
end
 
 
% get peristimulus times
%--------------------------------------------------------------------------
try
    
    % time window and bins for modelling
    %----------------------------------------------------------------------
    BSD.xY.Time = time(D, [], 'ms'); 
    T1          = BSD.options.Tdcm(1);
    T2          = BSD.options.Tdcm(2);
    [i, T1]     = min(abs(BSD.xY.Time - T1));
    [i, T2]     = min(abs(BSD.xY.Time - T2));
    
    % Time [ms] of down-sampled data
    %----------------------------------------------------------------------
    It          = [T1:DT:T2]';               % indices - bins
    BSD.xY.pst  = BSD.xY.Time(It);           % PST
    BSD.xY.It   = It;                        % Indices of time bins
    BSD.xY.dt   = DT/D.fsample;              % sampling in seconds
    Nb          = length(It);                % number of bins
    
catch
    errordlg('Please specify time window');
    error('')
end
 
% get frequency range
%--------------------------------------------------------------------------
try
    if length(BSD.options.Fdcm) > 2
        BSD.xY.Hz = BSD.options.Fdcm;
    else 
        Hz1     = BSD.options.Fdcm(1);          % lower frequency
        Hz2     = BSD.options.Fdcm(end);        % upper frequency
        BSD.xY.Hz  = fix(Hz1:Hz2);             % Frequencies
    end
catch
    pst     = BSD.xY.pst(end) - BSD.xY.pst(1);
    Hz1     = max(ceil(2*1000/pst),4);
    if Hz1 < 8
        Hz2 = 48;
    else
        Hz2 = 128;
    end
    BSD.xY.Hz  = fix(Hz1:Hz2);             % Frequencies
end


% Frequencies
%--------------------------------------------------------------------------
Nf         = length(BSD.xY.Hz);        % number of frequencies
Ne         = length(trial);            % number of trial types


% Cross spectral density for each trial type
%==========================================================================
condlabels = D.condlist;               % condition or trial type labels
BSD.xY.csd = cell(1,Ne);               % CSD for each condition

if 1/BSD.xY.dt < 2*BSD.xY.Hz (end)
    warning(['Requested frequency (%.1fHz) does not satisfy Nyquist ' ...
        'criterion for sampling rate (%.1fHz)\n'], BSD.xY.Hz (end),  1/BSD.xY.dt);
end

% w     = min(fix(2/BSD.xY.dt),Nb);      % window length (bins)
m     = 1;                             % retain principal mode
for i = 1:Ne
   
    % trial indices
    %----------------------------------------------------------------------
    c = D.indtrial(condlabels(trial(i)), 'GOOD');
    if ~BSD.M.noprint
        fprintf('\nevaluating CSD for condition %i\n',i)
    end
    
    % use only the first 512 trial
    %----------------------------------------------------------------------
    try c = c(1:512); end
    Nt    = length(c);
    
    % Get data
%     %----------------------------------------------------------------------
%     Nw    = max(8*(fix(Nb/w) - 1),1);
%     K     = zeros(Nf*Nm*Nm,Nw);
%     for k = 1:Nw
%         P     = zeros(Nf,Nm,Nm);
%         for j = 1:Nt
%             Iw  = It((1:w) + fix((k - 1)*w/8));
%             Y   = full(double(D(Ic,Iw,c(j))'*BSD.M.U));
%             mar = spm_mar(Y,32); % Note: Low order makes the spectrum smooth
%             mar = spm_mar_spectra(mar,BSD.xY.Hz,1/BSD.xY.dt);
%             P   = P + mar.P;
%         end
%         
%         % store
%         %------------------------------------------------------------------
%         K(:,k) = spm_vec(P/Nt);
%     end
%     
%     % retain principal eigenmode
%     %----------------------------------------------------------------------
%     [u,s,v]       = spm_svd(K,1);
%     P             = u(:,m)*s(m,m)*mean(v(:,m), 1)';

    P = 0;
    for j = 1:Nt
        Y   = full(double(D(Ic,:,c(j))'*BSD.M.U));
        
        P   = P + spm_csd(Y,BSD.xY.Hz,1/BSD.xY.dt);
    end
        
    BSD.xY.csd{i} = P; 
   
end

 
% place cross-spectral density in xY.y
%==========================================================================
BSD.xY.y    = BSD.xY.csd; 
BSD.xY.U    = BSD.M.U;
BSD.xY.code = condlabels(trial);
