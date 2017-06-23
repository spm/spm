function [images, outroot] = spm_eeg_convert2images(S)
% Convert M/EEG data to images for statistical analysis
% FORMAT [images, outroot] = spm_eeg_convert2images(S)
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
%                'source' (average over time and frequency)
%                'time x frequency' (average over channels)
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
%
% output:
%   images     - list of generated image files or objects
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak, James Kilner, Stefan Kiebel
% $Id: spm_eeg_convert2images.m 7125 2017-06-23 09:49:29Z guillaume $

SVNrev = '$Rev: 7125 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG conversion setup'); spm('Pointer','Watch');

if ~isfield(S, 'timewin'),      S.timewin  = [-Inf Inf];    end
if ~isfield(S, 'freqwin'),      S.freqwin  = [-Inf Inf];    end
if ~isfield(S, 'channels'),     S.channels = 'all';         end
if ~isfield(S, 'prefix'),       S.prefix = '';              end

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
    
    freqonset = D.frequencies(freqind(1));
    df = unique(diff(D.frequencies(freqind)));
    if length(df)> 1
        if (max(diff(df))/mean(df))>0.1
            df = unique(diff(log(D.frequencies(freqind))));
            if (max(diff(df))/mean(df))>0.1
                error('Irregular frequency spacing');
            else
                freqonset = log(D.frequencies(freqind(1)));
            end
        end
    end
    df = mean(df);
else
    df = 0;
end

isscalp = strncmpi('scalp', S.mode, 5);
ismesh  = false;

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
V     = [136 172]/n;        % voxel size

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
            0     0     df              freqonset;...
            0     0     0               1];
        N.mat(3,4) = N.mat(3,4) - N.mat(3,3);
        
        dat = file_array('', [n n length(freqind)], 'FLOAT32-LE');
        
    case 'scalp'
        avflag = [0 1 1];
        
        N.mat = [...
            V(1)  0     0               -C(1);...
            0     V(2)  0               -C(2);...
            0     0     1               1;...
            0     0     0               1];
        N.mat(3,4) = N.mat(3,4) - N.mat(3,3);
        
        dat = file_array('', [n n], 'FLOAT32-LE');
        
    case 'source'
        avflag = [0 1 1];
        
        g = D.sensors('SRC');
        
        if ~isempty(g)
            g = gifti(g);
            ismesh = true;
        end
        
        if isempty(g) || size(g.vertices, 1) ~= length(chanind)
            error('The number of source channels should match the number of mesh vertices.');
        end
        
    case 'time x frequency'
        if ~isTF
            error('This mode is only supported for TF datasets.');
        end
        
        avflag = [1 0 0];
        
        N.mat = [...
            df      0               0  freqonset;...
            0       1e3/D.fsample   0  time(D, timeind(1), 'ms');...
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
        N.mat(1, 4) = freqonset - N.mat(1, 1);
        
        dat = file_array('', [length(freqind) 1], 'FLOAT32-LE');
        
    case 'average'
        avflag = [1 1 1];
        
        N.mat(1, 4) = time(D, timeind(1), 'ms');
        
        if isTF
            N.mat(2, 4) = freqonset;
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

%-Create the output folder
%--------------------------------------------------------------------------
outroot    = [S.prefix spm_file(D.fname, 'basename')];
[sts, msg] = mkdir(D.path, outroot);
if ~sts,     error(msg); end
outroot     = fullfile(D.path, outroot);

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
    
    
    %-Make subdirectory for each condition
    %--------------------------------------------------------------------------
    condlabel = S.conditions{c};
    condlabel = condlabel(isstrprop(condlabel, 'alphanum') | ismember(condlabel, '_-'));
    if isempty(condlabel)
        condlabel = sprintf('condition%04d', c);
    end
    
    if ismesh
        res        = mkdir(outroot, ['condition_' condlabel]);
        
        save(g, fullfile(outroot, ['condition_' condlabel], [spm_file(D.fname, 'basename'), '.surf.gii']));
        
        G = gifti;
        G.private.metadata(1).name = 'SurfaceID';
        G.private.metadata(1).value =  fullfile(outroot, ['condition_' condlabel], [spm_file(D.fname, 'basename'), '.surf.gii']);
        
    else
        fname      = fullfile(outroot, ['condition_' condlabel '.nii']);
        
        cdat       = dat;
        cdat.fname = fname;
        cdat.dim   = [dat.dim ones(1, 3-length(dat.dim)) length(trialind)];
        N.dat      = cdat;
        create(N);
    end
    
    spm_progress_bar('Init', length(trialind),['Converting condition ',S.conditions{c}],'Trial');
    if length(trialind) > 100, Ibar = floor(linspace(1, length(trialind), 100)); else Ibar = 1:length(trialind); end
    
    for i = 1:length(trialind)
        
        Y = subsref(D, struct('type', '()', 'subs', {[dataind, {trialind(i)}]}));
        
        for m = sort(find(avflag), 1, 'descend')
            Y = mean(Y, m);
        end
        
        Y = squeeze(Y);
        
        if isscalp
            for j = 1:size(Y, 2)
                YY = NaN(n,n);
                YY(sub2ind([n n], x, y)) = griddata(Cel(:,1),Cel(:,2),...
                    double(Y(:, j)), x, y,'linear');
                
                switch length(N.dat.dim)
                    case 2
                        N.dat(:,:)     = YY;
                    case 3
                        N.dat(:,:, j)  = YY;
                    case 4
                        N.dat(:,:,j,i) = YY;
                    otherwise
                        error('Invalid output file');
                end
            end
            
            images{ind} = [fname ',' num2str(i)];
            
        elseif ismesh
            fname      = fullfile(outroot,  ['condition_' condlabel], ['condition_' condlabel '_' sprintf('_%04d',trialind(i)) '.gii']);
            
            G.cdata = Y(:);
            
            save(G, fname, 'ExternalFileBinary');
            
            images{ind} = fname;
        else
            if size(Y, 1) == 1
                Y = Y';
            end
            
            N.dat(:, :, 1, i) = Y;
            
            images{ind} = [fname ',' num2str(i)];
        end
        
        ind = ind + 1;
        
        if ismember(i, Ibar), spm_progress_bar('Set', i); end
    end
end


images = images(:);


%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
spm('FigName','M/EEG conversion: done'); spm('Pointer','Arrow');
