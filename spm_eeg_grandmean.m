function Do = spm_eeg_grandmean(S)
% average over multiple data sets
% FORMAT Do = spm_eeg_grandmean(S)
%
% S         - struct (optional)
% (optional) fields of S:
% D         - filenames (char matrix) of EEG mat-file containing epoched
%             data
% Dout      - filename (with or without path) of output file
% weighted  - average weighted by number of replications in inputs (1)
%             or not (0).
%
% Output:
% Do        - EEG data struct, result files are saved in the same
%                 directory as first input file.
%__________________________________________________________________________
%
% spm_eeg_grandmean averages data over multiple files. The data must have
% the same trialtype numbering and sampling rate. This function can be used
% for grand mean averaging, i.e. computing the average over multiple subjects.
% Missing event types and bad channels are taken into account properly.
% The output is written to a user-specified new file. The default name is
% the same name as the first selected input file, but prefixed with a 'g'.
% The output file is written to the current working directory.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_grandmean.m 3897 2010-05-21 15:06:50Z vladimir $

SVNrev = '$Rev: 3897 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG Grand Mean'); spm('Pointer','Watch');

%-Get MEEG object
%--------------------------------------------------------------------------
try
    S.D;
catch
    [S.D, sts] = spm_select([2 Inf], 'mat', 'Select M/EEG mat file');
    if ~sts, Do = []; return; end
end

%-Get parameters
%--------------------------------------------------------------------------
try
    S.Dout;
catch
    [filename, pathname] = uiputfile('*.mat', 'Select output file');
    S.Dout = fullfile(pathname, filename);
end

if ~isfield(S, 'weighted')
    S.weighted = spm_input('Weighted average?','+1','yes|no',[1 0], 1);
end

%-Load MEEG data
%--------------------------------------------------------------------------
D = cell(1, size(S.D, 1));
for i = 1:size(S.D, 1)
    try
        D{i} = spm_eeg_load(deblank(S.D(i, :)));
    catch
        error('Trouble reading files %s', deblank(S.D(i, :)));
    end
end
Nfiles = length(D);

%-Check dimension of the data files
%--------------------------------------------------------------------------
estr = [];

for i = 1:Nfiles
    flist = [];
    if ~strcmp(D{i}.type, 'evoked');
        flist = [flist ' ' D{i}.fname];
    end
    if ~isempty(flist)
       estr = ['The files' flist ' do not contain averaged (evoked) data.'];
    end
end

% This applies to # of channels, # of samples, # of conditions, sampling
% rate, and # of frequencies (if time-frequency data).
nc = zeros(Nfiles,1);
ns = zeros(Nfiles,1);
fs = zeros(Nfiles,1);
ne = zeros(Nfiles,1);
if strncmp(D{1}.transformtype, 'TF',2)
    nf = zeros(Nfiles,1);
end

for i = 1:Nfiles
    nc(i) = D{i}.nchannels;
    ns(i) = D{i}.nsamples;
    fs(i) = D{i}.fsample;
    ne(i) = D{i}.nconditions;
    if strncmp(D{1}.transformtype, 'TF',2)
        nf(i) = D{i}.nfrequencies;
    end
end

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
    %estr = [estr,'Data don''t have the same number of conditions.\n'];
end

if strncmp(D{1}.transformtype, 'TF',2)
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



%-Initialise output
%--------------------------------------------------------------------------
Do = D{1};

try
    [f1, f2, f3] = fileparts(S.Dout);
    Do.path = f1;
    fnamedat = [f2 '.dat'];
catch
    fnamedat = ['g' D{1}.fnamedat];
end

% how many different trial types and bad channels
types = {};
for i = 1:Nfiles
    types = unique([types, D{i}.condlist]);
end

% The order of the conditions will be consistent with the first file
[sel1, sel2] = spm_match_str(D{1}.condlist, types);
types = types([sel2, setdiff(1:length(types), sel2)]);

Ntypes = numel(types);

% how many repetitons per trial type
nrepl = zeros(Nfiles, Ntypes);
for i = 1:Nfiles
    for j = 1:numel(types)
        ind = D{i}.indtrial(types{j});
        if ~isempty(ind)
            nrepl(i, j) =  D{i}.repl(ind);
        end
    end 
end

if ~S.weighted
    nrepl = ones(Nfiles, Ntypes);
end

% generate new meeg object with new filenames
if strncmp(D{1}.transformtype, 'TF',2)
    Do = clone(Do, [spm_str_manip(S.Dout, 'r') '.dat'], [Do.nchannels Do.nfrequencies Do.nsamples Ntypes]);
else
    Do = clone(Do, [spm_str_manip(S.Dout, 'r') '.dat'], [Do.nchannels Do.nsamples Ntypes]);
end

% for determining bad channels of the grandmean
w = zeros(Do.nchannels, Ntypes);


%-Do the averaging
%--------------------------------------------------------------------------
spm_progress_bar('Init', Ntypes, 'responses averaged');
if Ntypes > 100, Ibar = floor(linspace(1, Ntypes, 100));
else Ibar = [1:Ntypes]; end

if strncmp(D{1}.transformtype, 'TF',2)
    for i = 1:Ntypes
        d = zeros(D{1}.nchannels, D{1}.nfrequencies, D{1}.nsamples);

        for j = 1:D{1}.nchannels

            for k = 1:Nfiles
                if ~ismember(j, D{k}.badchannels)
                    ind = D{k}.pickconditions(types{i});
                    if ~isempty(ind)
                        d(j, :, :) = d(j, :, :) + nrepl(k, i)*D{k}(j, :, :, ind);
                        w(j, i) = w(j, i) + nrepl(k, i);
                    end
                end
            end

            if w(j, i) > 0
                d(j, :, :) = d(j, :, :)/w(j, i);
            end

        end

        Do(1:Do.nchannels, 1:Do.nfrequencies, 1:Do.nsamples, i) = d;

        if ismember(i, Ibar), spm_progress_bar('Set', i); end

    end

else
    for i = 1:Ntypes
        d = zeros(D{1}.nchannels, D{1}.nsamples);

        for j = 1:D{1}.nchannels

            for k = 1:Nfiles
                if ~ismember(j, D{k}.badchannels)
                    ind = D{k}.pickconditions(types{i});
                    if ~isempty(ind)
                        d(j, :) = d(j, :) + nrepl(k, i)*D{k}(j, :, ind);
                        w(j, i) = w(j, i) + nrepl(k, i);
                    end
                end
            end

            if w(j, i) > 0
                d(j, :) = d(j, :)/w(j, i);
            end

        end

        Do(1:Do.nchannels, 1:Do.nsamples, i) = d;

        if ismember(i, Ibar), spm_progress_bar('Set', i); end

    end
end

spm_progress_bar('Clear');

%-Save Grand Mean to disk
%--------------------------------------------------------------------------
Do = type(Do, 'evoked');

bads = find(~any(w'));
if ~isempty(bads)
    Do = badchannels(Do, bads, 1);
else
    Do = badchannels(Do, 1:Do.nchannels, 0);
end

nrepl = sum(nrepl, 1);

Do = conditions(Do, [], types);
Do = repl(Do, [], nrepl);
Do = reject(Do, [], 0);
Do = trialonset(Do, [], []);
Do = Do.history('spm_eeg_grandmean', S, 'reset');

save(Do);

%-Display Grand Mean
%--------------------------------------------------------------------------
if ~isfield(S, 'review') || S.review
    spm_eeg_review(Do);
end

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG Grand Mean: done'); spm('Pointer','Arrow');
