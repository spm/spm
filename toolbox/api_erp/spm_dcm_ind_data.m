function DCM = spm_dcm_ind_data(DCM)
% gets time-frequency amplitude at specified sources for DCM
% FORMAT DCM = spm_dcm_ind_data(DCM)
% DCM    -  DCM structure
% requires
%
%    DCM.xY.Dfile
%    DCM.M.dipfit
%    DCM.options.Nmodes
%    DCM.options.trials
%    DCM.options.Tdcm
%    DCM.options.D    
%
% sets
%
%    DCM.xY.pst     - Peristimulus Time [ms] of time-frequency data
%    DCM.xY.dt      - sampling in seconds
%    DCM.xY.y       - concatenated induced response over sources
%    DCM.xY.xf      - induced response over sourcese
%    DCM.xY.It      - Indices of time bins
%    DCM.xY.Ic      - Indices of good channels
%    DCM.xY.Hz      - Frequency bins (for Wavelet transform)
%    DCM.xY.Mz      - Mean frequnecy response over trial and sources
%    DCM.xY.Rft     - wavelet coeficient
%    DCM.xY.Nm      - number of frequency modes
%    DCM.xY.U       - Frequnecy modes
%    DCM.xY.S       - and their singular values
%
%    DCM.options.h
%__________________________________________________________________________
% Stefan Kiebel, Karl friston
% $Id: spm_dcm_erp_data.m 668 2006-10-26 16:35:28Z karl $

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

% indices of EEG channel (excluding bad channels) and perstimulus times
%--------------------------------------------------------------------------
Ic              = setdiff(D.channels.eeg, D.channels.Bad);
Nc              = length(Ic);
DCM.M.dipfit.Ic = Ic;
DCM.xY.Ic       = Ic;
DCM.xY.Time     = 1000*[-D.events.start:D.events.stop]/D.Radc; % ms
DCM.xY.dt       = 1/D.Radc;

% options
%--------------------------------------------------------------------------
try
    Nm    = DCM.options.Nmodes;
catch
    errordlg('Please specify number of frequency modes');
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

% get peristimulus times
%--------------------------------------------------------------------------
try
    
    % time window and bins for modelling
    %----------------------------------------------------------------------
    T1      = DCM.options.Tdcm(1);
    T2      = DCM.options.Tdcm(2);
    [i, T1] = min(abs(DCM.xY.Time - T1));
    [i, T2] = min(abs(DCM.xY.Time - T2));
    
    % Time [ms] of downsampled data
    %----------------------------------------------------------------------
    DCM.xY.dt  = DT/D.Radc;                 % sampling in seconds
    B          = fix(0.064/DCM.xY.dt);      % 64ms boundary
    It         = [T1:DT:T2]';               % indices - bins
    Ns         = length(It);                % number of bins
    Nb         = length(It) - B - B;        % number of samples
    DCM.xY.pst = DCM.xY.Time(It(1 + B:end - B));        % PST
    DCM.xY.It  = It;                        % Indices of time bins

catch
    errordlg('Please specify time window');
    error('')
end

% get Morelet wavelets
%--------------------------------------------------------------------------
DCM.xY.Nm  = Nm;                       % number of frequency modes
DCM.xY.Hz  = 4:1:48;                   % Frequencies
ST         = DCM.xY.dt*1000;           % sampling interval
Nf         = length(DCM.xY.Hz);        % number of frequencies
Ne         = length(trial);            % number of ERPs
Nm         = DCM.xY.Nm;                % number of frequency modes
DCM.xY.Rft = linspace(5,5,Nf);         % wavelet coeficient


% get induced responses (use previous time-frequency results if possible) 
%==========================================================================
try
    if size(DCM.xY.xf,1) == Ne;
        if size(DCM.xY.xf,2) == Nm;
            if size(DCM.xY.xf{1},1) == Nb;
                DCM.xY.y  = spm_cat(DCM.xY.xf);
                return
            end
        end
    end
end

% high-pass filter
%--------------------------------------------------------------------------
h     = max(DCM.options.h,2);
T     = spm_dctmtx(Ns,h);
T     = speye(Ns,Ns) - T*T';

% create convolution matrices (Eucldian normalised with filtering)
%--------------------------------------------------------------------------
for i = 1:Nf
    W    = spm_eeg_morlet(DCM.xY.Rft(i), ST, DCM.xY.Hz(i));
    W    = convmtx(W{1}',Ns);
    N    = fix((size(W,1) - Nb)/2);
    W    = W([1:Nb] + N,:);
    M{i} = W*T;
end

% get gain matrix for source components
%==========================================================================

% parameterised lead field ECD given positions
%--------------------------------------------------------------------------
pos    = DCM.M.dipfit.L.pos;
mom    = [1  0  0;
          0  1  0;
          0  0  1];
Ng     = size(mom,2);                  % number of moments per source
Nr     = size(pos,2);                  % number of sources
G.Lmom = kron(ones(1,Nr),mom);
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
        c = find(D.events.code == D.events.types(i) & ~D.events.reject);
    else
        c = find(D.events.code == D.events.types(i));
    end
    Nt    = length(c);
    
    
    % Get data
    %----------------------------------------------------------------------
    Ny    = Nb*Nr*3;
    Y     = zeros(Ny,Nt);
    for j = 1:Nf
        f     = [1:Ny] + (j - 1)*Ny;
        for k = 1:Nt
            y      = abs(M{j}*D.data(Ic,It,c(k))'*MAP');
            Y(f,k) = y(:);
        end
        fprintf('\nevaluating %i Hz, condition %i',DCM.xY.Hz(j),i)
    end
    
    % weight with principal eigenvariate over trials
    %----------------------------------------------------------------------
    Y     = Y/normest(Y);
    u     = spm_svd(Y'*Y);
    u     = full(u(:,1)*sign(max(u(:,1))));
    Y     = reshape(Y*u,Nb,Nr*3,Nf);
    
    % sum time-frequency response over moments, normalise and remove mean
    %----------------------------------------------------------------------
    for j = 1:Nr
        Yk    = zeros(Nb,Nf);
        for k = 1:3
            Yk = Yk + squeeze(Y(:,j + k - 1,:))/Nt;
        end
        Mz{i,j} = mean(Yk);
        Yz{i,j} = Yk - ones(Nb,1)*Mz{i,j};
    end
end

% reduce to frequency modes
%==========================================================================

% find frequency modes (over time and sources)
%--------------------------------------------------------------------------
Y     = spm_cat(Yz(:));
[U S] = spm_svd(Y'*Y,0);
U     = U(:,1:Nm);

% project time-frequnecy data onto modes
%--------------------------------------------------------------------------
DCM.xY.xf = cell(Ne,Nm);
for i = 1:Ne
    for j = 1:Nr
        xf{j} = Yz{i,j}*U;
    end
    DCM.xY.xf(i,:) = spm_cell_swap(xf)';
end
DCM.xY.y  = spm_cat(DCM.xY.xf);
DCM.xY.U  = U;
DCM.xY.S  = S;
DCM.xY.Mz = Mz;



