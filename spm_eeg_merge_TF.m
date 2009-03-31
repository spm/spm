function Dout = spm_eeg_merge_TF(S)
% Concatenate epoched single trial files containing time-frequency data
% FORMAT D = spm_eeg_merge_TF(S)
%
% S         - optional input struct
% (optional) fields of S:
% D         - filename of EEG mat-file with continuous data
% recode    - a cell matrix of condition labels (rows: file, column: condition). 
%             
% Output:
% D         - EEG data struct (also written to files)
%
% concatenates epoched single trial files containing time-frequency data
%__________________________________________________________________________
%
% This function can be used to merge M/EEG files containing time-frequency 
% data to one file. This is useful whenever the data are distributed over 
% multiple files, but one wants to use all information in one file. For 
% example, when displaying data (SPM displays data from only one file at a 
% time), or merging information that has been measured in multiple sessions.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
% 
% Stefan Kiebel
% $Id: spm_eeg_merge_TF.m 3013 2009-03-31 13:52:57Z vladimir $

SVNrev = '$Rev: 3013 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG TF Merge'); spm('Pointer','Watch');

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select([2 inf], 'mat', 'Select M/EEG mat files');
    if ~sts, Dout = []; return; end
    S.D = D;
end


try
    for i = 1:size(D, 1)
        F{i} = spm_eeg_load(deblank(D(i, :)));
    end
    D = F;
catch
    error('Trouble reading files');
end

Nfiles = length(D);

if Nfiles < 2
    error('Need at least two files for merging');
end

Dout = D{1};

% Check input and determine number of new number of trial types
Ntrials = 0;
for i = 1:Nfiles
    if ~strncmp(D{i}.transformtype, 'TF', 2) % can be both TF or TFphase
        error('can only be used for time-frequency data');
    end
    
    % ascertain same number of channels, Nsamples and fsample, Nfrequencies
    if D{1}.nchannels ~= D{i}.nchannels
        error('Data don''t have the same number of channels.\nThere is a difference between files %s and %s.', D{1}.fname, D{i}.fname);
    end

    if D{1}.nsamples ~= D{i}.nsamples
        error('Data don''t have the same number of time points.\nThere is a difference between files %s and %s.', D{1}.fname, D{i}.fname);
    end

    if D{1}.fsample ~= D{i}.fsample
        error('Data don''t have the same sampling rate.\nThere is a difference between files %s and %s.', D{1}.fname, D{i}.fname);
    end

    if D{1}.frequencies ~= D{i}.frequencies
        error('Data don''t have the same frequencies.\nThere is a difference between files %s and %s.', D{1}.fname, D{i}.fname);
    end

    Ntrials = [Ntrials D{i}.ntrials];
end

% generate new meeg object with new filenames
[p, f, x] = fileparts(fnamedat(Dout));
Dout = clone(Dout, ['c' f x], [Dout.nchannels Dout.nfrequencies Dout.nsamples sum(Ntrials)]);
sDout = struct(Dout);

for i = 1:Nfiles

    clear Dtmp;
    Dtmp = D{i};
    cl = unique(Dtmp.conditions);

    % go to struct
    sDtmp = struct(Dtmp);

    % recode trial types
    try
        S.recode{i};
    catch
        for j = 1:Dtmp.nconditions
            S.recode{i}{j} = spm_input(sprintf('Labels: %s', spm_str_manip(Dtmp.fname, 'r')),...
                '+1', 's', sprintf('%s ', cl{j}));
        end
    end

    % recode labels
    code_new = {};
    for j = 1:Dtmp.nconditions
        [code_new{pickconditions(Dtmp, cl{j}, 0)}] = deal(S.recode{i}{j});
    end

    [sDtmp.trials.label] = deal(code_new{:});

    ind = [1:Ntrials(i+1)] + sum(Ntrials(1:i));

    [sDout.trials(ind).label] = deal(sDtmp.trials.label);
    
    % ignore updating onset timings
    try [sDout.trials(ind).onset] = deal(sDtmp.trials.onset); end
    try [sDout.trials(ind).repl] = deal(sDtmp.trials.repl); end

end

Dout = meeg(sDout);

% write files

% progress bar
spm_progress_bar('Init', Nfiles, 'Files merged');
if Nfiles > 100, Ibar = floor(linspace(1, Nfiles,100));
else Ibar = [1:Nfiles]; end

k = 0;
for i = 1:Nfiles

    ind = union(Dout.badchannels, D{i}.badchannels);
    if ~isempty(ind)
        Dout = badchannels(Dout, ind, 1);
    end

    % write trial-wise to avoid memory mapping error
    for j = 1:D{i}.ntrials
        k = k + 1;
        Dout(1:Dout.nchannels, 1:Dout.nfrequencies, 1:Dout.nsamples, k) =  D{i}(:,:,:,j);
        Dout = reject(Dout, k, reject(D{i}, j));
    end

    if ismember(i, Ibar), spm_progress_bar('Set', i); end


end

Dout = Dout.history('spm_eeg_merge_TF', S);
save(Dout);

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
spm('FigName','M/EEG TF merge: done'); spm('Pointer','Arrow');
