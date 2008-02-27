function DCM = spm_dcm_ind_data(DCM)
% gets time-frequency amplitude at specified sources for DCM
% FORMAT DCM = spm_dcm_ind_data(DCM)
% DCM    -  DCM structure
% requires
%
%    DCM.xY.Dfile
%    DCM.Lpos
%    DCM.options.Nmodes
%    DCM.options.trials
%    DCM.options.Tdcm
%    DCM.options.Fdcm
%    DCM.options.D
%    DCM.options.Rft
%    DCM.options.h    
%
% sets
%
%    DCM.xY.pst     - Peristimulus Time [ms] of time-frequency data
%    DCM.xY.dt      - sampling in seconds
%    DCM.xY.y       - concatenated induced response over sources
%    DCM.xY.xf      - induced response over sources
%    DCM.xY.It      - Indices of time bins
%    DCM.xY.Ic      - Indices of good channels
%    DCM.xY.Hz      - Frequency bins (for Wavelet transform)
%    DCM.xY.Mz      - Mean frequnecy response over trial and sources
%    DCM.xY.Rft     - wavelet coeficient
%    DCM.xY.Nm      - number of frequency modes
%    DCM.xY.U       - Frequnecy modes
%    DCM.xY.S       - and their singular values
%
%    DCM.xY.y{i}(k,l)    = l-th region X frequency mode (fast over regions)
%                          k-th time-bin
%                          i-th trial
%
%    DCM.xY.xf{i,j}(k,l) = l-th frequnecy mode
%                          k-th time-bin
%                          j-th region
%                          i-th trial
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_ind_data.m 1174 2008-02-27 20:22:30Z karl $

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
D = spm_eeg_ldata(Dfile);

% indices of EEG channel (excluding bad channels)
%--------------------------------------------------------------------------
Ic              = setdiff(D.channels.eeg, D.channels.Bad);
Nc              = length(Ic);
DCM.xY.Ic       = Ic;

% options
%--------------------------------------------------------------------------
try
    Nm    = DCM.options.Nmodes;
catch
    errordlg('Please specify number of frequency modes');
    error('')
end
try
    h     = DCM.options.h;
catch
    errordlg('Please number of DCT components');
    error('')
end
try
    DT    = DCM.options.D;
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
    DCM.xY.Rft = DCM.options.Rft;                                        
catch
    % default wavelet number
    %----------------------------------------------------------------------
    DCM.xY.Rft = 5;
end

% check data are not oversampled (< 4ms)
%--------------------------------------------------------------------------
if DT/D.Radc < 0.004
    DT            = ceil(0.004*D.Radc);
    DCM.options.D = DT;
end


% get peristimulus times
%--------------------------------------------------------------------------
try
    
    % time window and bins for modelling
    %----------------------------------------------------------------------
    DCM.xY.Time = 1000*[-D.events.start:D.events.stop]/D.Radc; % ms
    T1          = DCM.options.Tdcm(1);
    T2          = DCM.options.Tdcm(2);
    [i, T1]     = min(abs(DCM.xY.Time - T1));
    [i, T2]     = min(abs(DCM.xY.Time - T2));
    
    % Time [ms] of downsampled data
    %----------------------------------------------------------------------
    Ns          = length(DCM.xY.Time);       % number of bins
    It          = [T1:DT:T2]';               % indices - bins
    Is          = [1:Ns]';                   % indices - samples
    DCM.xY.pst  = DCM.xY.Time(It);           % PST
    DCM.xY.It   = It;                        % Indices of time bins
    DCM.xY.dt   = DT/D.Radc;                 % sampling in seconds
    Nb          = length(It);                % number of bins
    
catch
    errordlg('Please specify time window');
    error('')
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

if (Hz2 - Hz1) > 64, HzD = 2; else, HzD = 1; end

% get Morelet wavelets
%--------------------------------------------------------------------------
DCM.xY.Hz  = Hz1:HzD:Hz2;              % Frequencies
DCM.xY.Nm  = Nm;                       % number of frequency modes
dt         = 1000/D.Radc;              % sampling interval (ms)
Nf         = length(DCM.xY.Hz);        % number of frequencies
Nr         = size(DCM.C,1);            % number of sources
Ne         = length(trial);            % number of ERPs
Nm         = DCM.xY.Nm;                % number of frequency modes

% get induced responses (use previous time-frequency results if possible) 
%==========================================================================
try
    if size(DCM.xY.xf,1) == Ne;
        if size(DCM.xY.xf,2) == Nr;
            if size(DCM.xY.xf{1},1) == Nb;
                if size(DCM.xY.xf{1},2) == Nm;
                    if size(DCM.xY.U,1) == length(DCM.xY.Hz)
                        DCM.xY.y = spm_cond_units(spm_cat(spm_cell_swap(DCM.xY.xf),2));
                        return
                    end
                end
            end
        end
    end
end

% high-pass filter (detrend)
%--------------------------------------------------------------------------
T     = spm_orthpoly(Ns,h);
T     = speye(Ns,Ns) - T*T';

% create convolution matrices (Eucldian normalised with filtering)
%--------------------------------------------------------------------------
for i = 1:Nf
    W    = spm_eeg_morlet(DCM.xY.Rft, dt, DCM.xY.Hz(i));
    N    = fix(length(W{1})/2);
    W    = convmtx(W{1}',Ns);
    W    = W(It + N,:);
    M{i} = W*T;
end

% get gain matrix for source components
%==========================================================================

% parameterised lead field ECD given positions
%--------------------------------------------------------------------------
try
    pos = DCM.Lpos;
catch
    pos = DCM.M.dipfit.L.pos;
end
mom    = [1  0  0;
          0  1  0;
          0  0  1];
Ng     = size(mom,2);                  % number of moments per source
Nr     = size(pos,2);                  % number of sources
G.L    = kron(ones(1,Nr),mom);
G.Lpos = kron(pos,ones(1,Ng));
L      = spm_erp_L(G,DCM.M);
MAP    = pinv(L);

% add (spatial filtering) re-referencing to MAP projector
%--------------------------------------------------------------------------
R     = speye(Nc,Nc) - ones(Nc,1)*pinv(ones(Nc,1));
MAP   = MAP*R;


% Wavelet amplitudes for each (projected) source
%==========================================================================
for i = 1:Ne;
    
    % trial indices
    %----------------------------------------------------------------------
    if isfield(D.events,'reject')
        c = find(D.events.code == D.events.types(trial(i)) & ~D.events.reject);
    else
        c = find(D.events.code == D.events.types(trial(i)));
    end
    
    % use only the first 512 trial
    %----------------------------------------------------------------------
    try c = c(1:512); end
    Nt    = length(c);
    
    
    % Get data
    %----------------------------------------------------------------------
    Ny    = Nb*Nr*3;
    Y     = zeros(Ny,Nt);
    for j = 1:Nf
        f     = [1:Ny] + (j - 1)*Ny;
        for k = 1:Nt
            y      = abs(M{j}*D.data(Ic,Is,c(k))'*MAP');
            Y(f,k) = y(:);
        end
        fprintf('\nevaluating %i Hz, condition %i (%i trials)',DCM.xY.Hz(j),i,Nt)
    end
    
    % weight with principal eigenvariate over trials (c.f., averaging)
    %----------------------------------------------------------------------
    u     = spm_svd(Y'*Y);
    u     = full(u(:,1)*sign(max(u(:,1))));
    Y     = reshape(Y*u,Nb,Nr*3,Nf);
    
    % sum time-frequency response over moments and remove baseline
    %----------------------------------------------------------------------
    for j = 1:Nr
        Yk    = zeros(Nb,Nf);
        for k = 1:3
            Yk = Yk + squeeze(Y(:,j + k - 1,:))/Nt;
        end
        Mz{i,j} = Yk(1,:);
        Yz{i,j} = Yk - ones(Nb,1)*Mz{i,j};
    end
end

% reduce to frequency modes
%==========================================================================

% find frequency modes (over time and sources)
%--------------------------------------------------------------------------
Y          = spm_cat(Yz(:));
[U S]      = spm_svd(Y'*Y,0);
U          = U(:,1:Nm);

% project time-frequency data onto modes
%--------------------------------------------------------------------------
DCM.xY.xf = cell(Ne,Nr);
for i = 1:Ne
    for j = 1:Nr
        DCM.xY.xf{i,j} = Yz{i,j}*U;
    end
end
DCM.xY.y    = spm_cond_units(spm_cat(spm_cell_swap(DCM.xY.xf),2));
DCM.xY.U    = U;
DCM.xY.S    = S;
DCM.xY.code = D.events.code(trial);



