function D = spm_eeg_weight_epochs(S);
% computes contrasts over trials or trial types.
% FORMAT D = spm_eeg_weight_epochs(S)
%
% S         - optional input struct
% (optional) fields of S:
% D         - filename of EEG mat-file with epoched data
% c         - contrast matrix, each row computes a contrast of the data
%
% Output:
% D         - EEG data struct (also written to files)
%_______________________________________________________________________
%
% spm_eeg_weight_epochs computes contrasts of data, over epochs of data. The
% input is a single MEEG file.
% The argument c must have dimensions Ncontrasts X Nepochs, where Ncontrasts is
% the number of contrasts and Nepochs the number of epochs, i.e. each row of c
% contains one contrast vector. The output
% is a MEEG file with Ncontrasts epochs. The typical use is to compute,
% for display purposes, contrasts like the difference or interaction
% between trial types in channel space.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel, Rik Henson
% $Id: spm_eeg_weight_epochs.m 1237 2008-03-21 14:54:07Z stefan $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG averaging setup',0);

try
    D = S.D;
catch
    D = spm_select(1, '.*\.mat$', 'Select EEG mat file');
end

P = spm_str_manip(D, 'H');

try
    D = spm_eeg_load(D);
catch
    error('Trouble reading file %s', D);
end

try
    c = S.c;
catch
    c = spm_input('Enter contrasts', '1', 'x', '', inf, eye(D.ntrials));
end

if ~isempty(D.repl)
    try
        WeightAve = S.WeightAve;
    catch
        WeightAve = spm_input('Weight by num replications?', '+1', 'yes|no', [1 0]);
    end
else
    WeightAve = 0;
end

% here should be a sanity check of c


spm('Pointer', 'Watch'); drawnow;

Ncontrasts = size(c, 1);

% generate new meeg object with new filenames
Dnew = newdata(D, ['m' fnamedat(D)], [D.nchannels D.nsamples Ncontrasts], dtype(D));

spm_progress_bar('Init', Ncontrasts, 'Contrasts computed'); drawnow;
if Ncontrasts > 100, Ibar = floor(linspace(1, Ncontrasts, 100));
else Ibar = [1:Ncontrasts]; end

for i = 1:Ncontrasts

    if WeightAve
        p = find(c(i,:) == 1);
        if ~isempty(p)
            r = D.repl(p);
            c(i,p) = r/sum(r);
        end

        p = find(c(i,:) == -1);
        if ~isempty(p)
            r = D.repl(p);
            c(i,p) = -r/sum(r);
        end
    end

    disp(['Contrast ', mat2str(i),': ', mat2str(c(i,:),3)])

    d = zeros(D.nchannels, D.nsamples);

    for j = 1:D.nchannels
        d(j, :) = c(i,:) * squeeze(D(j, :, :))';
    end

    Dnew = putdata(Dnew, 1:Dnew.nchannels, 1:Dnew.nsamples, i, d);

    newrepl(i) = sum(D.repl(find(c(i,:)~=0)));

    if ismember(i, Ibar)
        spm_progress_bar('Set', i);
        drawnow;
    end


end

spm_progress_bar('Clear');

% jump outside methods to reorganise trial structure
sD = struct(Dnew);
sD.trials = sD.trials(1:Ncontrasts);

for i = 1:Ncontrasts
    sD.trials(i).code = sprintf('contrast %d', i);
    try sD.trials(i) = rmfield(sD.trials(i), 'onset'); end
    try sD.trials(i) = rmfield(sD.trials(i), 'reject'); end
    try sD.trials(i) = rmfield(sD.trials(i), 'blinks'); end
end

try [sD.trials.repl] = deal(newrepl); end

Dnew = meeg(sD);

save(Dnew);

spm('Pointer', 'Arrow');
