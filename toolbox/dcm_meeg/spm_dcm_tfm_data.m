function DCM = spm_dcm_tfm_data(DCM)
% gets cross-spectral density data-features using a wavelet transform
% FORMAT DCM = spm_dcm_tfm_data(DCM)
% DCM    -  DCM structure
% requires
%
%    DCM.xY.Dfile        - name of data file
%    DCM.M.U             - channel subspace
%    DCM.options.trials  - trial to evaluate
%    DCM.options.Tdcm    - time limits
%    DCM.options.Fdcm    - frequency limits
%    DCM.options.D       - Down-sampling
%
% sets
%
%    DCM.xY.pst     - Peristimulus Time [ms] sampled
%    DCM.xY.dt      - sampling in seconds [s] (down-sampled)
%    DCM.xY.U       - channel subspace
%    DCM.xY.y       - cross spectral density over channels
%    DCM.xY.csd     - cross spectral density over channels
%    DCM.xY.erp     - event-related average over channels
%    DCM.xY.It      - Indices of time bins
%    DCM.xY.Ic      - Indices of good channels
%    DCM.xY.Hz      - Frequency bins
%    DCM.xY.code    - trial codes evaluated
%    DCM.xY.scale   - scalefactor applied to data
%    DCM.xY.Rft     - Wavelet number or ratio of frequency to time
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_tfm_data.m 6234 2014-10-12 09:59:10Z karl $
 
% Set defaults and Get D filename
%-------------------------------------------------------------------------
try
    Dfile = DCM.xY.Dfile;
catch
    errordlg('Please specify data and trials');
    error('')
end
 
% load D
%--------------------------------------------------------------------------
try
    D = spm_eeg_load(Dfile);
catch
    try
        [dummy,f]    = fileparts(Dfile);
        D            = spm_eeg_load(f);
        DCM.xY.Dfile = fullfile(pwd,f);
    catch
        try
            [f,p]        = uigetfile('*.mat','please select data file');
            name         = fullfile(p,f);
            D            = spm_eeg_load(name);
            DCM.xY.Dfile = fullfile(name);
        catch
            warndlg([Dfile ' could not be found'])
            return
        end
    end
end

% options
%--------------------------------------------------------------------------
try, DT    = DCM.options.D;      catch, DT    = 1;             end
try, trial = DCM.options.trials; catch, trial = D.nconditions; end
try, Rft   = DCM.options.Rft;    catch, Rft   = 8;             end
try, han   = DCM.options.han;    catch, han   = 0;             end
try, h     = DCM.options.h;      catch, h     = 1;             end
 

% Modality 
%--------------------------------------------------------------------------
if ~isfield(DCM.xY, 'modality')
    [mod, list] = modality(D, 0, 1);
    
    if isequal(mod, 'Multimodal')
        qstr = 'Only one modality can be modelled at a time. Please select.';
        if numel(list) < 4
            
            % Nice looking dialog
            %--------------------------------------------------------------
            options         = [];
            options.Default = list{1};
            options.Interpreter = 'none';
            DCM.xY.modality = questdlg(qstr, 'Select modality', list{:}, options);
            
        else
            
            % Ugly but can accommodate more buttons
            %--------------------------------------------------------------
            ind = menu(qstr, list);
            DCM.xY.modality = list{ind};
        end
    else
        DCM.xY.modality = mod;
    end
end
 
% channel indices (excluding bad channels)
%--------------------------------------------------------------------------
DCM.xY.Ic = D.indchantype(DCM.xY.modality,'GOOD');
Ic        = DCM.xY.Ic;

% ensure spatial modes have been computed (see spm_dcm_csd)
%-------------------------------------------------------------------------
if ~isfield(DCM.M,'U')
    DCM.M.U = spm_dcm_eeg_channelmodes(DCM.M.dipfit,DCM.options.Nmodes);
end
if size(DCM.M.U,1) ~= length(Ic);
    DCM.M.U = spm_dcm_eeg_channelmodes(DCM.M.dipfit,DCM.options.Nmodes);
end
Nm          = size(DCM.M.U,2);

 
% check data are not oversampled (< 4ms)
%--------------------------------------------------------------------------
if DT/D.fsample < 0.004
    DT            = ceil(0.004*D.fsample);
    DCM.options.D = DT;
end


% get frequency range
%--------------------------------------------------------------------------
try
    Hz1     = DCM.options.Fdcm(1);          % lower frequency
    Hz2     = DCM.options.Fdcm(2);          % upper frequency
catch
    pst     = DCM.xY.pst(end) - DCM.xY.pst(1);
    Hz1     = max(ceil(2*1000/pst),4);
    if Hz1 < 8;
        Hz2 = 48;
    else
        Hz2 = 128;
    end
end
 
 
% Frequencies
%--------------------------------------------------------------------------
Hz  = fix(Hz1:Hz2);                         % Frequencies
Nf  = length(Hz);                           % number of frequencies
Ne  = length(trial);                        % number of ERPs

 
% get peristimulus times
%--------------------------------------------------------------------------
try
    % padding for time frequency analysis
    %----------------------------------------------------------------------
    Np          = ceil(Rft/min(Hz)*D.fsample/DT/2);
    
    % time window and bins for modelling
    %----------------------------------------------------------------------
    DCM.xY.Time = time(D, [], 'ms');
    T1          = DCM.options.Tdcm(1);
    T2          = DCM.options.Tdcm(2);
    [dummy, T1] = min(abs(DCM.xY.Time - T1));
    [dummy, T2] = min(abs(DCM.xY.Time - T2));
    P1          = T1 - Np*DT;
    P2          = T2 + Np*DT;
    
    % Time [ms] of down-sampled data
    %----------------------------------------------------------------------
    It          = (T1:DT:T2)';               % indices - bins
    Ip          = (P1:DT:P2)';               % indices - bins
    DCM.xY.pst  = DCM.xY.Time(It);           % PST
    DCM.xY.It   = It;                        % indices of time bins
    DCM.xY.dt   = DT/D.fsample;              % sampling in seconds
    Nb          = length(It);                % number of bins
    Ns          = length(Ip);                % number of padded samples
    Ib          = (1:Nb) + Np;               % indices of padded samples
    
catch
    errordlg('Please specify time window');
    error('')
end

% confounds - DCT:
%--------------------------------------------------------------------------
if h == 0
    X0 = sparse(Ns,1);
else
    X0 = spm_dctmtx(Ns,h);
end
R      = speye(Ns) - X0*X0';
 
% hanning (omit second residualization for very long time-series)
%--------------------------------------------------------------------------
if han
    if Ns < 2048
        R = R*sparse(diag(hanning(Ns)))*R;
    else
        R = sparse(diag(hanning(Ns)))*R;
    end
end
 
 
% Cross spectral density for each trial type
%==========================================================================
cond  = D.condlist;
erp   = cell(1,Ne);
csd   = cell(1,Ne);
for e = 1:Ne;
    
    % trial indices
    %----------------------------------------------------------------------
    c = D.indtrial(cond(trial(e)), 'GOOD');
    
    % use only the first 512 trial
    %----------------------------------------------------------------------
    try c = c(1:512); end
    
    
    % evoked response
    %----------------------------------------------------------------------
    Nt    = length(c);
    Y     = zeros(Nb + 2*Np,Nm,Nt);
    P     = zeros(Nb,Nf,Nm,Nm);
    Q     = zeros(Nb,Nf,Nm,Nm);
    
    for k = 1:Nt
        Y(:,:,k) = R*D(Ic,Ip,c(k))'*DCM.M.U;
    end
 
    
    % ERP or average
    %----------------------------------------------------------------------
    A = mean(Y,3);
    
    % induced response
    %----------------------------------------------------------------------
    for k = 1:Nt
        
        fprintf('\nevaluating condition %i (trial %i)',e,k)
        G     = spm_morlet(Y(:,:,k) - A,Hz*DCM.xY.dt,Rft);
        for i = 1:Nm
            for j = 1:Nm
                P(:,:,i,j) = (G(Ib,:,i).*conj(G(Ib,:,j)));
            end
        end
        Q = Q + P;
    end
    
    % normalise induced responses in relation variance about ERP
    %--------------------------------------------------------------------------
    Vm    = mean(mean(squeeze(var(Y,[],3))));
    Vs    = mean(diag(squeeze(mean(squeeze(mean(Q))))));
    Q     = Vm*Q/Vs;
    
    % store
    %----------------------------------------------------------------------
    csd{e} = Q;
    erp{e} = A(Ib,:);
    
end
 
% place cross-spectral density in xY.y
%==========================================================================
try
    scale = DCM.xY.scale;
catch
    scale = 1;
end

DCM.xY.erp   = spm_unvec(spm_vec(erp)*(scale^1),erp);
DCM.xY.csd   = spm_unvec(spm_vec(csd)*(scale^2),csd);
DCM.xY.y     = DCM.xY.csd;
DCM.xY.U     = DCM.M.U;
DCM.xY.scale = scale;
DCM.xY.Hz    = Hz;
DCM.xY.Rft   = Rft;
 
return
 
% plot responses
%--------------------------------------------------------------------------
spm_dcm_tfm_response(DCM.xY,DCM.xY.pst,DCM.xY.Hz);
