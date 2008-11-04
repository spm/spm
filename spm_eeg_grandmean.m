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
% $Id: spm_eeg_grandmean.m 2436 2008-11-04 10:46:27Z stefan $

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


%% Check dimension of the data files
% This applies to # of channels, # of samples, # of conditions, sampling
% rate, and # of frequencies (if time-frequency data).
nc = zeros(Nfiles,1);
ns = zeros(Nfiles,1);
fs = zeros(Nfiles,1);
ne = zeros(Nfiles,1);
if strcmp(D{1}.transformtype, 'TF')
    nf = zeros(Nfiles,1);
end

for i = 1:Nfiles
    nc(i) = D{i}.nchannels;
    ns(i) = D{i}.nsamples;
    fs(i) = D{i}.fsample;
    ne(i) = D{i}.nconditions;
    if strcmp(D{1}.transformtype, 'TF')
        nf(i) = D{i}.nfrequencies;
    end
end

estr = [];
unc = unique(nc);
if length(unc)~=1
    ind = zeros(Nfiles,1);
    fna = cell(length(unc),1);
    for i=1:Nfiles
        ind(i) = find(unc==nc(i));
        fna{ind(i)} = [fna{ind(i)},'-',D{i}.fname,'\n'];
    end
    fprintf('\n')
    for i=1:length(unc)
        fprintf(['\nThose files have ',num2str(unc(i)),' channels:\n',...
            fna{i}])
    end
    estr = [estr,'Data don''t have the same number of channels.\n'];
end

uns = unique(ns);
if length(uns)~=1
    ind = zeros(Nfiles,1);
    fna = cell(length(uns),1);
    for i=1:Nfiles
        ind(i) = find(uns==ns(i));
        fna{ind(i)} = [fna{ind(i)},'-',D{i}.fname,'\n'];
    end
    fprintf('\n')
    for i=1:length(uns)
        fprintf(['\nThose files have ',num2str(uns(i)),' time points:\n',...
            fna{i}])
    end
    estr = [estr,'Data don''t have the same number of time points.\n'];
end

ufs = unique(fs);
if length(ufs)~=1
    ind = zeros(Nfiles,1);
    fna = cell(length(ufs),1);
    for i=1:Nfiles
        ind(i) = find(ufs==fs(i));
        fna{ind(i)} = [fna{ind(i)},'-',D{i}.fname,'\n'];
    end
    fprintf('\n')
    for i=1:length(ufs)
        fprintf(['\nThose files have a sampling rate of ',num2str(ufs(i)),' Hz:\n',...
            fna{i}])
    end
    estr = [estr,'Data don''t have the same sampling rate.\n'];
end

une = unique(ne);
if length(une)~=1
    ind = zeros(Nfiles,1);
    fna = cell(length(une),1);
    for i=1:Nfiles
        ind(i) = find(une==ne(i));
        fna{ind(i)} = [fna{ind(i)},'-',D{i}.fname,'\n'];
    end
    fprintf('\n')
    for i=1:length(unc)
        fprintf(['\nThose files have ',num2str(une(i)),' conditions:\n',...
            fna{i}])
    end
    estr = [estr,'Data don''t have the same number of conditions.\n'];
end

if strcmp(D{1}.transformtype, 'TF')
    unf = unique(nf);
    if length(unf)~=1
        ind = zeros(Nfiles,1);
        fna = cell(length(unf),1);
        for i=1:Nfiles
            ind(i) = find(unf==nf(i));
            fna{ind(i)} = [fna{ind(i)},'-',D{i}.fname,'\n'];
        end
        fprintf('\n')
        for i=1:length(unf)
            fprintf(['\nThose files have ',num2str(unf(i)),' frequency bins:\n',...
                fna{i}])
        end
        estr = [estr,'Data don''t have the same number of frequency bins.\n'];
    end
end

% send message error (if any)
if ~isempty(estr)
    error(estr)
    return
else
    fprintf('Ok: All data files have the same dimensions.\n')
end



%% Initialization
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
    Do = clone(Do, [spm_str_manip(S.Pout, 'r') '.dat'], [Do.nchannels Do.nfrequencies Do.nsamples Ntypes]);
else
    Do = clone(Do, [spm_str_manip(S.Pout, 'r') '.dat'], [Do.nchannels Do.nsamples Ntypes]);
end

% for determining bad channels of the grandmean
w = zeros(Do.nchannels, Ntypes);


%% Do the averaging
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

bads = find(~any(w'));
if ~isempty(bads)
    Do = badchannels(Do, bads, 1);
else
    Do = badchannels(Do, 1:Do.nchannels, 0);
end

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

Do = meeg(sD);
Do = Do.history('spm_eeg_grandmean', S);

save(Do);

if ~isfield(S, 'review') || S.review
    spm_eeg_review(Do);
end

spm('Pointer', 'Arrow');