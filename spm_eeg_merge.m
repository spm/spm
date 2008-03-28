function Dout = spm_eeg_merge(S);
% FORMAT D = spm_eeg_merge(S)
%
% S         - optional input struct
% (optional) fields of S:
% D         - filename of EEG mat-file with continuous data
% recode    - a cell vector with one cell for each file. A vector in each cell
%             codes the new event type. This allows to either recode events
%             or keep their old codes.
% Output:
% D         - EEG data struct (also written to files)
%
% concatenates epoched single trial files
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_merge.m 1262 2008-03-28 09:28:42Z stefan $

% Changed to allow recoding of first file (though obviously makes some of
% loop redundant!)          Doris Eckstein
% Corrected checks for reject, repl and Bad fields         RH 8/1/08

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG merge',0);

try
    D = S.D;
catch
    D = spm_select([2 inf], '\.mat$', 'Select M/EEG mat files');
end

P = spm_str_manip(D(1,:), 'H');

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

spm('Pointer', 'Watch');

Dout = D{1};

% Check input and determine number of new number of trial types
Ntrials = [];
for i = 1:Nfiles
    % ascertain same number of channels, Nsamples and fsample
    if D{1}.nchannels ~= D{i}.nchannels
        error('Data don''t have the same number of channels.\nThere is a difference between files %s and %s.', D{1}.fname, D{i}.fname);
    end

    if D{1}.nsamples ~= D{i}.nsamples
        error('Data don''t have the same number of time points.\nThere is a difference between files %s and %s.', D{1}.fname, D{i}.fname);
    end

    if D{1}.fsample ~= D{i}.fsample
        error('Data don''t have the same sampling rate.\nThere is a difference between files %s and %s.', D{1}.fname, D{i}.fname);
    end

    Ntrials = [Ntrials D{i}.ntrials];
end

% generate new meeg object with new filenames
Dout = newdata(Dout, ['c' fnamedat(Dout)], [Dout.nchannels Dout.nsamples sum(Ntrials)], dtype(Dout));


for i = 1:Nfiles

    clear Dtmp;
    Dtmp = D{i};

    % go to struct
    sDout = struct(Dout);
    sDtmp = struct(Dtmp);

    % recode trial types
    try
        S.recode{i};
    catch
        S.recode{i} = spm_input(sprintf('Labels: %s', spm_str_manip(Dtmp.fname, 'r')),...
            '+1', 's+', sprintf('%s ', Dtmp.conditionlabels), Dtmp.nconditions);
    end

    % recode labels
    cl = cellstr(Dtmp.conditionlabels);
    code_new = {};
    for j = 1:Dtmp.nconditions
        [code_new(pickconditions(Dtmp, cl{j}))] = deal(S.recode{1});
    end

    [sDtmp.trials.label] = deal(code_new{:});

    if i > 1
        ind = [1:Ntrials(i)] + sum(Ntrials(i-1));

        [sDout.trials(ind).label] = deal(sDtmp.trials.label);
        % update onset timings: only relevant for continuous data?
        try [sDout.trials(ind).onset] = deal(sDtmp.trials.onset + Dtmp.fsample/1000*Dtmp.nsamples); end
        try [sDout.trials(ind).repl] = deal(sDtmp.trials.repl); end
    end
end

Dout = meeg(sDout);

% write files

% progress bar
spm_progress_bar('Init', Nfiles, 'Files merged'); drawnow;
if Nfiles > 100, Ibar = floor(linspace(1, Nfiles,100));
else Ibar = [1:Nfiles]; end

k = 0;
for i = 1:Nfiles

    Dout = putbadchannels(Dout, intersect(Dout.badchannels, D{i}.badchannels), 1);

    % write trial-wise to avoid memory mapping error
    for j = 1:D{i}.ntrials
        k = k + 1;
        Dout = putdata(Dout, 1:Dout.nchannels, 1:Dout.nsamples, k, D{i}(:,:,j));
    end

    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end


end

save(Dout);
spm_progress_bar('Clear');
spm('Pointer', 'Arrow');



