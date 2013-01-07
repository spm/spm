function images = spm_eeg_convert2images(S)
% Convert M/EEG data to images for statistical analysis
% FORMAT images = spm_eeg_convert2images(S)
%
% S                   - input structure (optional)
%  fields of S:
%   D          - MEEG object or filename of M/EEG mat-file with
%                epoched data
%
%   mode       - type of images to generate one of:
%                'scalp x time'
%                'scalp x frequency' (average over time)
%                'scalp' (average over time and frequency)
%                'time x frequency' (avrage over channels)
%                'time' (1D average over channels, frequency)
%                'frequency' (1D average over channels, time)
%                'average' (average over all dimensions to get a single
%                           number)
%
%   conditions - cell array of condition labels (default: convert all
%                conditions)
%   timewin    - time window to retain (in PST ms)
%   freqwin    - frequency window to retain (for TF datasets)
%   channels   - cell array of channel labels, modality or 'all'.
%
%   prefix     - prefix for the folder containing the images (default: none)
%   virtual    - if 0 (default) images will always be written to disk
%                1  return virtual handles if possible
%
% output:
%   images     - list of generated image files or objects
%__________________________________________________________________________
% Copyright (C) 2005-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, James Kilner, Stefan Kiebel
% $Id: spm_eeg_convert2images.m 5177 2013-01-07 11:36:08Z vladimir $

SVNrev = '$Rev: 5177 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG conversion setup'); spm('Pointer','Watch');

if ~isfield(S, 'timewin'),      S.timewin  = [-Inf Inf];    end
if ~isfield(S, 'freqwin'),      S.freqwin  = [-Inf Inf];    end
if ~isfield(S, 'channels'),     S.channels = 'all';         end
if ~isfield(S, 'virtual'),      S.virtual = 0;              end

D = spm_eeg_load(S.D);

if ~isfield(S, 'conditions') || isempty(S.conditions),  S.conditions = D.condlist;  end
if ~iscell(S.conditions), S.conditions = {S.conditions};                            end

if strcmp(D.type, 'continuous')
    error('Continuous data are not supported');
end

isTF = strncmpi(D.transformtype,'TF',2);

timeind = D.indsample(1e-3*(min(S.timewin))):D.indsample(1e-3*(max(S.timewin)));
if isempty(timeind) || any(isnan(timeind))
    error('Selected time window is invalid.');
end

if isTF
    freqind = D.indfrequency(min(S.freqwin)):D.indfrequency(max(S.freqwin));
    if isempty(freqind) || any(isnan(freqind))
        error('Selected frequency window is invalid.');
    end
    
    df = unique(diff(D.frequencies(freqind)));
    if length(df)> 1;
        error('Irregular frequency spacing');
    end
else
    df = 0;
end

isscalp = strncmpi('scalp', S.mode, 5);

chanind = setdiff(D.selectchannels(S.channels), D.badchannels);

if isempty(chanind)
    error('No channels were selected');
end

modality = unique(D.chantype(chanind));
if length(modality)>1
    error('All channels should be of the same type. Process each modality separately.');
end

if isequal(modality, 'MEGPLANAR') && isscalp
    error('Planar channels should be combined before creating scalp maps');
end


N     = nifti;
N.mat = eye(4);
N.mat_intent = 'Aligned';

C     = [68  100];  % origin
n     = 32;         % dimension (make the default from SPM8 constant here)
V     = C/n;        % voxel size

switch S.mode
    case 'scalp x time'
        avflag = [0 1 0];
        
        N.mat = [...
            V(1)  0     0               -C(1);...
            0     V(2)  0               -C(2);...
            0     0     1e3/D.fsample   time(D, timeind(1), 'ms');...
            0     0     0               1];
        N.mat(3,4) = N.mat(3,4) - N.mat(3,3);
        
        dat = file_array('', [n n length(timeind)], 'FLOAT32-LE');
        
    case 'scalp x frequency'
        if ~isTF
            error('This mode is only supported for TF datasets.');
        end
        
        avflag = [0 0 1];
        
        N.mat = [...
            V(1)  0     0               -C(1);...
            0     V(2)  0               -C(2);...
            0     0     df              D.frequencies(freqind(1));...
            0     0     0               1];
        N.mat(3,4) = N.mat(3,4) - N.mat(3,3);
        
        dat = file_array('', [n n length(freqind)], 'FLOAT32-LE');
        
    case 'scalp'
        avflag = [0 1 1];
        
        N.mat = [...
            V(1)  0     0               -C(1);...
            0     V(2)  0               -C(2);...
            0     0     1               1e3*mean(D.time(timeind));...
            0     0     0               1];
        N.mat(3,4) = N.mat(3,4) - N.mat(3,3);
        
        dat = file_array('', [n n 1], 'FLOAT32-LE');
        
    case 'time x frequency'
        if ~isTF
            error('This mode is only supported for TF datasets.');
        end
        
        avflag = [1 0 0];
        
        N.mat = [...
            df      0               0  D.frequencies(freqind(1));...
            0       1e3/D.fsample   0  time(D, 1, 'ms');...
            0       0               1  0;...
            0       0               0  1];
        N.mat(1,4) = N.mat(1,4) - N.mat(1,1);
        N.mat(2,4) = N.mat(2,4) - N.mat(2,2);
        
        dat = file_array('', [length(freqind) length(timeind)], 'FLOAT32-LE');
        
    case 'time'
        avflag = [1 1 0];
        
        N.mat(1, 1) = 1e3/D.fsample;
        N.mat(1, 4) = time(D, timeind(1), 'ms') - N.mat(1, 1);
        
        dat = file_array('', [length(timeind) 1], 'FLOAT32-LE');
        
    case 'frequency'
        if ~isTF
            error('This mode is only supported for TF datasets.');
        end
        
        avflag = [1 0 1];
        
        N.mat(1, 1) = df;
        N.mat(1, 4) = D.frequencies(freqind(1)) - N.mat(1, 1);
        
        dat = file_array('', [length(freqind) 1], 'FLOAT32-LE');
        
    case 'average'
        avflag = [1 1 1];
        
        N.mat(1, 4) = time(D, timeind(1), 'ms');
        
        if isTF
            N.mat(2, 4) = D.frequencies(freqind(1));
        end
        
        dat = file_array('', [1 1], 'FLOAT32-LE');
    otherwise
        error('Unsupported mode.');
end

avflag = logical(avflag);

if isTF
    dataind = {chanind, freqind, timeind};
else
    avflag = avflag([1 3]);
    dataind = {chanind, timeind};
end

% Force saving when the virtual mode is impossible
virtual = S.virtual && ~(isscalp || any(cellfun(@length, dataind(avflag)) > 1));

if virtual ~= S.virtual
    warning('Images must be written to disk for this case. Overriding the user setting');
end

if ~virtual
    outroot    = [S.prefix spm_file(D.fname, 'basename')];
    [sts, msg] = mkdir(D.path, outroot);
    if ~sts,     error(msg); end
    outroot     = fullfile(D.path, outroot);
end

if isscalp
    [Cel, x, y] = spm_eeg_locate_channels(D, n, chanind);
end

images = {};
ind    = 1;
for c = 1:numel(S.conditions)
    trialind = D.indtrial(S.conditions{c}, 'GOOD');
    
    if isempty(trialind)
        warning(['No good trials found for condition ' S.conditions{c}]);
        continue;
    end          
    
    if ~virtual
        %-Make subdirectory for each condition
        %--------------------------------------------------------------------------
        outdir = ['condition_' S.conditions{c}];
        [sts, msg] = mkdir(outroot, outdir);
        if ~sts, error(msg); end
        outdir = fullfile(outroot, outdir);
    end
    
    spm_progress_bar('Init', length(trialind),['Converting condition ',S.conditions{c}],'Trial');
    if length(trialind) > 100, Ibar = floor(linspace(1, length(trialind), 100)); else Ibar = 1:length(trialind); end
    
    for i = 1:length(trialind)
        if ~virtual
            if strcmp(D.type, 'single')
                % single trial data
                fname = [sprintf('trial%04d', trialind(i)) spm_file_ext];
            else
                % evoked data
                fname = ['average' spm_file_ext];
            end
            
            fname = fullfile(outdir, fname);
            
            images{ind} = fname;
            ind = ind + 1;
            
            dat.fname = fname;
            N.dat = dat;
            create(N);
            
            Y = subsref(D, struct('type', '()', 'subs', {[dataind, {trialind(i)}]}));
            
            for m = sort(find(avflag), 1, 'descend')
                Y = mean(Y, m);
            end
            
            Y = squeeze(Y);
            
            if isscalp
                for j = 1 : size(Y, 2)
                    YY = NaN(n,n);
                    YY(sub2ind([n n], x, y)) = griddata(Cel(:,1),Cel(:,2),...
                        double(Y(:, j)), x, y,'linear');
                    N.dat(:,:,j) = YY;
                end                
            else                               
                if size(Y, 1) == 1
                    Y = Y';
                end
                
                N.dat(:, :) = Y;
            end
        end
        
        if ismember(i, Ibar), spm_progress_bar('Set', i); end
    end
end


if ~virtual
    images = char(images);
end

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
spm('FigName','M/EEG conversion: done'); spm('Pointer','Arrow');
