function Do = spm_eeg_grandmean(S)
% average over multiple data sets
% FORMAT Do = spm_eeg_grandmean(S)
%
% S         - struct (optional)
% (optional) fields of S:
% P         - filenames (char matrix) of EEG mat-file containing epoched
%             data
% Pout      - filename (with or without path) of output file
%
% Output:
% Do        - EEG data struct, result files are saved in the same
%                 directory as first input file.
%_______________________________________________________________________
%
% spm_eeg_grandmean averages data over multiple files. The data must have
% the same trialtype numbering and sampling rate. This function can be used for
% grand mean averaging, i.e. computing the average over multiple subjects.
% Missing event types and bad channels are taken into account properly.
% The output is written to a user-specified new file. The default name is
% the same name as the first selected input file, but prefixed with a 'g'.
% The output file is written to the current working directory.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_grandmean.m 2073 2008-09-10 10:25:39Z vladimir $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG grandmean setup', 0);

try
    P = S.P;
catch
    P = spm_select(inf, '\.mat$', 'Select EEG mat files');
    S.P = P;
end

try
    S.Pout;
catch
    [filename, pathname] = uiputfile('*.mat', 'Select output file');
    S.Pout = fullfile(pathname, filename);
end

D = cell(1, size(P, 1));
for i = 1:size(P, 1)
    try
        D{i} = spm_eeg_load(deblank(P(i, :)));
    catch
        error('Trouble reading files %s', deblank(P(i, :)));
    end
end

Nfiles = length(D);

% check input
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

    if strcmp(D{1}.transformtype, 'TF')
        if D{1}.nfrequencies ~= D{i}.nfrequencies
            error('Data don''t have the same number of frequencies.\nThere is a difference between files %s and %s.', D{1}.fname, D{i}.fname);
        end

    end
end

% output
Do = D{1};

try
    [f1, f2, f3] = fileparts(S.Pout);
    Do.path = f1;
    fnamedat = [f2 '.dat'];
catch
    fnamedat = ['g' D{1}.fnamedat];
end

spm('Pointer', 'Watch'); drawnow;

% how many different trial types and bad channels
types = {};
for i = 1:Nfiles
    types = unique([types, D{i}.conditions]);
end

Ntypes = numel(types);

% how many repetitons per trial type
repl = zeros(1, Ntypes);
for i = 1:Nfiles
    cl{i} = D{i}.conditions;
    for j = 1:D{i}.nconditions
        ind = strmatch(cl{i}{j}, types);
        repl(ind) =  repl(ind) + D{i}.repl(j);
    end
end

% generate new meeg object with new filenames
if strcmp(D{1}.transformtype, 'TF')
    Do = clone(Do, [spm_str_manip(S.Pout, 'tr') '.dat'], [Do.nchannels Do.nfrequencies Do.nsamples Ntypes]);
else
    Do = clone(Do, [spm_str_manip(S.Pout, 'tr') '.dat'], [Do.nchannels Do.nsamples Ntypes]);
end

% for determining bad channels of the grandmean
w = zeros(Do.nchannels, Ntypes);

spm_progress_bar('Init', Ntypes, 'responses averaged'); drawnow;
if Ntypes > 100, Ibar = floor(linspace(1, Ntypes, 100));
else Ibar = [1:Ntypes]; end

if strcmp(D{1}.transformtype, 'TF')
    for i = 1:Ntypes
        d = zeros(D{1}.nchannels, D{1}.nfrequencies, D{1}.nsamples);

        for j = 1:D{1}.nchannels

            for k = 1:Nfiles
                if ~ismember(j, D{k}.badchannels)
                    ind = strmatch(types(i), cl{k}(:));
                    if ~isempty(ind)
                        d(j, :, :) = d(j, :, :) + D{k}(j, :, :, ind);
                        w(j, i) = w(j, i) + 1;
                    end
                end
            end

            if w(j, i) > 0
                d(j, :, :) = d(j, :, :)/w(j, i);
            end

        end

        Do(1:Do.nchannels, 1:Do.nfrequencies, 1:Do.nsamples, i) = d;

        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end

    end

else
    for i = 1:Ntypes
        d = zeros(D{1}.nchannels, D{1}.nsamples);

        for j = 1:D{1}.nchannels

            for k = 1:Nfiles
                if ~ismember(j, D{k}.badchannels)
                    ind = strmatch(types(i), cl{k}(:));
                    if ~isempty(ind)
                        d(j, :) = d(j, :) + D{k}(j, :, ind);
                        w(j, i) = w(j, i) + 1;
                    end
                end
            end

            if w(j, i) > 0
                d(j, :) = d(j, :)/w(j, i);
            end

        end

        Do(1:Do.nchannels, 1:Do.nsamples, i) = d;

        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end

    end
end
spm_progress_bar('Clear');

Do = type(Do, 'grandmean');

badchannels(Do, find(~any(w')), 1);
% jump to struct to make a few changes
sD = struct(Do);
for i = 1:Ntypes
    sD.trials(i).code = types{i};
    sD.trials(i).repl = repl(i);
end
% [sD.trials.repl] = deal(repl);
try sD.trials = rmfield(sD.trials, 'reject'); end
try sD.trials = rmfield(sD.trials, 'onset'); end
try sD.trials = rmfield(sD.trials, 'blinks'); end
sD.trials = sD.trials(1:Ntypes);
D = meeg(sD);

D = D.history('spm_eeg_grandmean', S);

save(D);

if ~isfield(S, 'review') || S.review
    spm_eeg_review(D);
end

spm('Pointer', 'Arrow');