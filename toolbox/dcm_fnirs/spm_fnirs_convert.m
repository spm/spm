function [D] = spm_fnirs_convert(S)
% Convert fNIRS data formats to SPM format, to perform DCM-fNIRS analysis 
% FORMAT [D] = spm_fnirs_convert(S) 
%
% S.fnirs  - file name of optical density changes 
% S.ch     - file name of source and detector position 
% S.fgreen - file name of Green's function
%
% D        - structure array of fNIRS data, which is an input of
%            spm_dcm_fnirs_specify.m 
%--------------------------------------------------------------------------
% D.y: optical density changes [# samples x (# of channels x # of waves)]
% Various fNIRS data formats are converted to matrix of optical density
% changes using the SPM-fNIRS toolbox (web) 
%
% D.pos: MNI coordinates of optical source and detectors 
% Source and detector positions in subject spaces are transformed to MNI
% coordinates using NFRI functions (web) 
% 
% D.fgreen: file name of Green's function. 
% Green's function can be estimated using MMC, MCX, or Toast software. 
% Estimated Green's function of the photon fluence should be saved as
% specific formats to be used in DCM-fNIRS analysis 
% (See spm_fnirs_sensitivity.m). 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny & Sungho Tak
% $Id: spm_fnirs_convert.m 6422 2015-04-23 16:55:52Z spm $

spm_input('Read fNIRS data files: ... ', 1, 'd');

%-Specify file names, and load data 
%--------------------------------------------------------------------------

%- Optical density changes 
%--------------------------------------------------------------------------
if ~nargin || ~isfield(S, 'fnirs') 
    [S.fnirs, sts] = spm_select([1 Inf], '^*.*\.mat$', 'Select fNIRS data files.');
    if ~sts, D = []; return; end
    
    nfiles = size(S.fnirs, 1); % number of files 
    nch = []; % number of channels 
    nwav = []; % number of wavelength 
    y = []; % optical density changes 
    fs = []; % sampling rate
    
    for i = 1:nfiles
        % if multiple data are selected, time courses are concatenated. 
        [fpath, fname, fext] = fileparts(S.fnirs(i,:));
        str = sprintf('fNIRS file%i: %s%s', i, fname, fext); 
        spm_input(str, '+1', 'd');
        
        load(S.fnirs(i,:)); 
        ns = size(nirs_data.OD, 1); 
        y{i,1} = reshape(nirs_data.OD, ns, []); 
        nch = [nch nirs_data.nch]; 
        nwav = [nwav size(nirs_data.OD, 3)]; 
        fs = [fs nirs_data.fs]; 
    end
    
    % check if formats of fNIRS data are consistent. 
    error_flag = sum(diff(nch)) + sum(diff(nwav)) + sum(diff(round(fs))); 
    if error_flag ~= 0, error('Formats of fNIRS data are not consistent.'); 
    else, nch = nch(1); nwav = nwav(1); fs = fs(1); end;
end

%- Source and detector positions 
%--------------------------------------------------------------------------
if ~nargin || ~isfield(S, 'ch') 
    [S.ch, sts] = spm_select(1, '^*.*\.mat$', 'Select optode position file.'); 
    if ~sts, D = []; return; end 
    [fpath, fname, fext] = fileparts(S.ch); 
    spm_input(sprintf('Optode position: %s%s', fname, fext), '+1', 'd'); 
    
    try
        load(S.ch); % optical probe positions
    catch
        error('Cannot load source and detector position file.');
    end
end

%- Green's function 
%--------------------------------------------------------------------------
if ~nargin || ~isfield(S, 'fgreen')
    [S.fgreen, sts] = spm_select(1, 'mat', 'Select a file of Greens function.'); 
    if ~sts, D = []; return; end
    spm_input(sprintf('Greens function: %s', S.fgreen), '+1', 'd');
end

%- Preprocess fNIRS measurements (if necessary) 
%--------------------------------------------------------------------------
spm_input('Preprocess time series of fNIRS data: ...', 1, 'd');

% Temporal filtering 
ftype = spm_input('Temporal filtering?', '+1', 'Butterworth IIR|No'); 
if strcmpi(ftype, 'Butterworth IIR') 
    addpath(fullfile(spm('Dir'), 'external', 'fieldtrip', 'preproc'));

    % cutoff frequencies 
    % [0 0.008]: very low-frequency drifts 
    % [0.12 0.35]: respiration 
    % [0.7 1.5]: cardiac pulsation 
    fcutoff = spm_input('cut-off frequencies [lower upper]:', '+1', 'r', '[0 0.008; 0.12 0.35; 0.7 1.5]');
    
    forder = 5; 
    ftype = 'but';
    nf = size(fcutoff, 1); 
    fdir = 'twopass';
    
    for i = 1:nfiles
        for j = 1:nf
            if fcutoff(j, 1) == 0 % high pass filter
                y{i,1} = ft_preproc_highpassfilter(y{i,1}',fs,fcutoff(j, 2), forder, ftype, fdir)';
            elseif fcutoff(j, 2) == Inf % low pass filter
                y{i,1} = ft_preproc_lowpassfilter(y{i,1}',fs,fcutoff(j, 1), forder, ftype, fdir)';
            else % band-stop filter
                y{i,1} = ft_preproc_bandstopfilter(y{i,1}',fs,fcutoff(j, :), forder, ftype, fdir)';
            end
        end
    end
end
    
% Concatenate time series 
str = [];
if nfiles > 1,
    fons = zeros(1, nfiles); 
    
    str{1,1} = 'fNIRS time series were concatenated:';
    fons = 0; 
    for i = 1:nfiles, 
        str = [str; sprintf('Onset of data %i (s): %4.1f', i, fons)];
        fons = (fons + size(y{i},1))./fs; 
    end
end

y = cell2mat(y);
for i = 1:size(str, 1), spm_input(str{i,1}, '+1', 'd'); end

% Temporal averaging 
ftype = spm_input('Temporal averaging?', '+1', 'y/n'); 
if strcmpi(ftype, 'y') 
    nw = spm_input('number of time window', 2, 'w1'); 
    wy = []; 
    for i = 1:nw
        str = sprintf('onsets of time window %i (s):', i); 
        ons = spm_input(str, '+1', 'r', ' ', [Inf 1]);
        mdur = mean(diff(ons));
        str = sprintf('window length (s):'); 
        dur = spm_input(str, '+1', 'r', mdur, 1); 
        
        ons = round(ons.*fs) + 1; dur = round(dur .* fs); 
        wy{i,1} = spm_fnirs_wavg(y, ons, dur);
    end
    y = cell2mat(wy); 
end

y = reshape(y, [], nch, nwav); 

%- Save data as structure array D 
%--------------------------------------------------------------------------
D.y = y; 
D.dt = 1./fs;
D.eps = nirs_data.eps;
D.nch = nch;
D.nwav = nwav; 

D.pos.s = pos.s; % source positions (MNI coordinates) 
D.pos.d = pos.d; % detector positions
D.pos.ch_sd = pos.ch_sd; % channel configuration  

D.fgreen = S.fgreen; 
