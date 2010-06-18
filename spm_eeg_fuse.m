function D = spm_eeg_fuse(S)
% Fuse MEG and EEG datasets to create a multimodal dataset
% FORMAT D = spm_eeg_fuse(S)
%
% S           - input structure (optional)
% (optional) fields of S:
%   S.D       - character array containing filenames of M/EEG mat-files
%
% D        - MEEG object (also written to disk, with a 'u' prefix)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
%
% Vladimir Litvak
% $Id: spm_eeg_fuse.m 3937 2010-06-18 14:41:43Z vladimir $

SVNrev = '$Rev: 3937 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG fuse'); spm('Pointer','Watch');

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select([2 Inf], 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

%-Load MEEG data
%--------------------------------------------------------------------------
F = cell(1,size(D,1));
try
    for i = 1:size(D, 1)
        F{i} = spm_eeg_load(deblank(D(i, :)));
    end
    D = F;
catch
    error('Trouble reading files');
end

Nfiles = numel(D);

if Nfiles < 2
    error('Need at least two files for fusion.');
end

%-Check input and determine number of new number of trial types
%--------------------------------------------------------------------------
channels = {};
modalities = {};
isTF   =  strncmpi(D{1}.transformtype,'TF',2); % TF and TFphase

for i = 1:Nfiles
    if ~isequal(D{i}.transformtype, D{1}.transformtype)
        error(['The datasets do not contain the same kind of data.\n'...
            'There is a difference between files\n\t%s\nand\n\t%s.'], ...
            D{1}.fname, D{i}.fname);
    end
    
    if D{1}.ntrials ~= D{i}.ntrials
        error(['Data don''t have the same number of trials.\n' ...
            'There is a difference between files\n\t%s\nand\n\t%s.'], ...
            D{1}.fname, D{i}.fname);
    end
    
    if ~isequal(D{1}.conditions, D{i}.conditions)
        error(['Data don''t have the same condition labels.\n' ...
            'There is a difference between files\n\t%s\nand\n\t%s.'], ...
            D{1}.fname, D{i}.fname);
    end
    
    if D{1}.nsamples ~= D{i}.nsamples
        error(['Data don''t have the same number of time points.\n' ...
            'There is a difference between files\n\t%s\nand\n\t%s.'], ...
            D{1}.fname, D{i}.fname);
    end
    
    if D{1}.fsample ~= D{i}.fsample
        error(['Data don''t have the same sampling rate.\n' ...
            'There is a difference between files\n\t%s\nand\n\t%s.'], ...
            D{1}.fname, D{i}.fname);
    end
    
    if isTF &&  ~isequal(D{1}.frequencies, D{i}.frequencies)
        error(['Data don''t have the same frequencies.\n' ...
            'There is a difference between files\n\t%s\nand\n\t%s.'], ...
            D{1}.fname, D{i}.fname);
    end
    
    if ~isempty(intersect(channels, D{i}.chanlabels))
        error(['Files to be fused should not have overlapping channel sets.\n' ...
            'There is an overlap in channel sets between files\n\t%s\nand\n\t%s.'], ...
            D{1}.fname, D{i}.fname);
    else
        channels = [channels, D{i}.chanlabels];
    end
    
    modalities = [modalities {D{i}.modality}];
end

Nchannels = numel(channels);

%-Arrange modalities in a fixed order
%--------------------------------------------------------------------------
[sel1, sel2] = spm_match_str({'Multimodal', 'MEG', 'EEG', 'LFP', 'Other'}, modalities);
D = D(sel2);
modalities = modalities(sel2);

%-Generate new meeg object with new filenames
%--------------------------------------------------------------------------
Dout = D{1};
[p, f, x] = fileparts(fnamedat(Dout));

if ~isTF
    Dout = clone(Dout, fullfile(pwd, ['u' f x]), [Nchannels Dout.nsamples Dout.ntrials]);
else
    Dout = clone(Dout, fullfile(pwd, ['u' f x]), [Nchannels Dout.nfrequencies Dout.nsamples Dout.ntrials]);
end


%-Write files
%--------------------------------------------------------------------------

spm_progress_bar('Init', Dout.ntrials, 'Trials written');
if Dout.ntrials > 100, Ibar = floor(linspace(1, Dout.ntrials, 100));
else Ibar = [1:Dout.ntrials]; end

for t = 1:Dout.ntrials
    k = 1;
    for i = 1:Nfiles
        if ~isTF
            Dout(k:(k+D{i}.nchannels-1), :, t) =  D{i}(:, :, t);
        else
            Dout(k:(k+D{i}.nchannels-1), :, :, t) =  D{i}(:, :, :, t);
        end
    k = k + D{i}.nchannels;
    end
    if ismember(t, Ibar), spm_progress_bar('Set', t); end
end

spm_progress_bar('Clear');

% Update the header data separately to only do it once
k = 1;
rejected = Dout.reject;
for i = 1:Nfiles
    rejected = rejected | D{i}.reject;
    
    Dout = chanlabels(Dout, k:(k+D{i}.nchannels-1), chanlabels(D{i}));
    Dout = chantype(Dout, k:(k+D{i}.nchannels-1), chantype(D{i}));
    if ~isempty(badchannels(D{i}))
        Dout = badchannels(Dout, k+badchannels(D{i})-1, 1);
    end
    Dout = units(Dout, k:(k+D{i}.nchannels-1), units(D{i}));
    Dout = coor2D(Dout, k:(k+D{i}.nchannels-1), coor2D(D{i}));
    
    k = k + D{i}.nchannels;
end

if any(rejected)
    Dout = reject(Dout, [], rejected);
end

%-Set sensor locations
%--------------------------------------------------------------------------
[junk, newmodalities] = modality(Dout, 1);
if ismember('MEG', newmodalities) && isempty(Dout.sensors('MEG'))
    for i = 1:numel(modalities)
        if isequal(modalities{i}, 'MEG') && ~isempty(D{i}.sensors('MEG'))
            Dout = sensors(Dout, 'MEG', D{i}.sensors('MEG'));
            Dout = fiducials(Dout, D{i}.fiducials);
            break;
        end
    end
end

if ismember('EEG', newmodalities) && isempty(Dout.sensors('EEG'))
    warning('Assigning default EEG sensor locations. Reload individual locations if necessary.');
    S1 = [];
    S1.D = Dout;
    S1.task = 'defaulteegsens';
    S1.updatehistory = 0;
    Dout = spm_eeg_prep(S1);
end

%-Remove previous inversions.
%--------------------------------------------------------------------------
if isfield(Dout, 'inv')
    Dout = rmfield(Dout, 'inv');
end

%-Save new M/EEG data
%--------------------------------------------------------------------------
Dout = Dout.history(mfilename, S, 'reset');
save(Dout);

D = Dout;

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
spm('FigName','M/EEG fuse: done'); spm('Pointer','Arrow');
