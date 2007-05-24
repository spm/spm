function DCM = spm_dcm_ind_data(DCM)
% gets time-frequency amplitude at specified sources for DCM
% FORMAT DCM = spm_dcm_ind_data(DCM)
% DCM    -  DCM structure
% requires
%
%    DCM.xY.Dfile
%    DCM.M.dipfit
%    DCM.options.trials
%    DCM.options.Tdcm
%    DCM.options.D    
%
% sets
%
%    DCM.xY.Time    - Time [ms] of downsampled data
%    DCM.xY.dt      - sampling in seconds
%    DCM.xY.y       - concatenated induced response over sources
%    DCM.xY.xf      - induced response over sourcese
%    DCM.xY.It      - Indices of time bins
%    DCM.xY.Ic      - Indices of good channels
%    DCM.xY.Hz      - Frequency bins (for Wavelet transform)
%    DCM.xY.Mz      - Mean frequnecy response over trial and sources
%    DCM.xY.Rft     - wavelet coeficient
%    DCM.xY.Nm      - number of requency modes
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

% if MEG, store grad struct in D.channels
%--------------------------------------------------------------------------
try
    DCM.xY.grad = D.channels.grad;
end

% get induced responses (use previous time-frequency results if possible) 
%==========================================================================
try
     DCM.xY.y  = spm_cat(DCM.xY.xf);
     return
end
    
% get Morelet wavelets
%--------------------------------------------------------------------------
DCM.options.h = 3;
DCM.xY.Nm     = 4;                        % number of frequnecy modes
DCM.xY.Rft    = 6;                        % wavelet coeficient
DCM.xY.Hz     = 4:1:64;                   % Frequencies
ST            = DCM.xY.dt*1000;           % sampling interval
Nf            = length(DCM.xY.Hz);

% high-pass filter
%--------------------------------------------------------------------------
h     = DCM.options.h;
T     = spm_dctmtx(Ns,h);
T     = speye(Ns,Ns) - T*T';

% create convolution matrices with reflecting boundaries (and filtering)
%--------------------------------------------------------------------------
M     = spm_eeg_morlet(DCM.xY.Rft, ST, DCM.xY.Hz, 32);
for i = 1:Nf
    M{i} = convmtx(M{i}',Ns);
    N    = fix((size(M{i},1) - Ns)/2);
    W    = M{i}([1:Ns] + N,:)*T;
    M{i} = diag(1./sqrt(sum(abs(W).^2,2)))*W;
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
for i = 1:length(trial);
    
    % trial indices
    %----------------------------------------------------------------------
    if isfield(D.events,'reject')
        c = find(D.events.code == D.events.types(i) & ~D.events.reject);
    else
        c = find(D.events.code == D.events.types(i));
    end
    Nt    = length(c);

    % Cycle over frequencies
    %----------------------------------------------------------------------
    Y     = zeros(Nf,Ns,Nr*3);
    for k = 1:Nf
        Yk    = zeros(Ns,Nr*3);
        for j = 1:Nt
            y  = MAP*D.data(Ic,It,c(j));
            Yk = Yk + abs(M{k}*y');
        end
        Y(k, :, :) = Yk;
        fprintf('\nevaluating frequency %i, condition %i',k,i)
    end
    
    % time-frequency repsones for trials and sources (summing over moments)
    %----------------------------------------------------------------------
    for j = 1:Nr
        Yk    = zeros(Ns,Nf);
        for k = 1:3
            Yk = Yk + Y(:,:,j + k - 1)'/Nt;
        end
        Mz{i,j} = mean(Yk);
        Yz{i,j} = spm_detrend(Yk);
    end
end


% reduce to frequency modes
%==========================================================================

% find frequency modes (over time and sources)
%--------------------------------------------------------------------------
Y     = spm_cat(Yz(:));
[U S] = spm_svd(Y'*Y,0);
U     = U(:,1:DCM.xY.Nm);

% project time-frequnecy data onto modes
%--------------------------------------------------------------------------
for i = 1:length(trial)
    for j = 1:Nr
        xf{j} = Yz{i,j}*U;
    end
    DCM.xY.xf(i,:) = spm_cell_swap(xf)';
end
DCM.xY.y  = spm_cat(DCM.xY.xf);
DCM.xY.U  = U;
DCM.xY.S  = S;
DCM.xY.Mz = Mz;



