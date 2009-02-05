function D = spm_eeg_weight_epochs_TF(S)
% Compute contrasts over trials or trial types, for time-frequency data
% FORMAT D = spm_eeg_weight_epochs(S)
%
% S         - optional input struct
% (optional) fields of S:
% D         - filename of EEG mat-file with epoched data
% c         - contrast matrix, each row computes a contrast of the data
% WeightAve - flag whether average should be weighted by number of
%             replications (yes (1), no (0))
% Output:
% D         - EEG data struct (also written to files)
%__________________________________________________________________________
%
% spm_eeg_weight_epochs computes contrasts of data, over epochs of data, 
% for time-frequency data. The input is a single MEEG file.
% The argument c must have dimensions Ncontrasts X Nepochs, where Ncontrasts is
% the number of contrasts and Nepochs the number of epochs, i.e. each row of c
% contains one contrast vector. The output
% is a M/EEG file with Ncontrasts epochs. The typical use is to compute,
% for display purposes, contrasts like the difference or interaction
% between trial types in channel space. Another possible use is remove
% trials from the data file, by using a contrast that contains zeros for
% the to be removed file.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_weight_epochs_TF.m 2696 2009-02-05 20:29:48Z guillaume $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG averaging setup',0);

try
    D = S.D;
catch
    D = spm_select(1, '.*\.mat$', 'Select EEG mat file');
    S.D = D;
end

P = spm_str_manip(D, 'H');

try
    D = spm_eeg_load(D);
catch
    error('Trouble reading file %s', D);
end

if ~strncmp(D.transformtype, 'TF', 2)
    error('Use this function for time-frequency data only');
end

try
    c = S.c;
catch
    c = spm_input('Enter contrasts', '1', 'x', '', inf, eye(D.ntrials));
    S.c = c;
end

if ~isempty(D.repl)
    try
        WeightAve = S.WeightAve;
    catch
        WeightAve = spm_input('Weight by num replications?', '+1', 'yes|no', [1 0]);
        S.WeightAve = WeightAve;
    end
else
    WeightAve = 0;
end

% here should be a sanity check of c


spm('Pointer', 'Watch'); drawnow;

Ncontrasts = size(c, 1);

% generate new meeg object with new filenames
Dnew = clone(D, ['m' fnamedat(D)], [D.nchannels D.nfrequencies D.nsamples Ncontrasts]);

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

    d = zeros(D.nchannels, D.nfrequencies, D.nsamples);

    for j = 1:D.nchannels
        for f = 1:D.nfrequencies
            d(j, f, :) = c(i,:) * squeeze(D(j, f, :, :))';
        end
    end

    Dnew(1:Dnew.nchannels, 1:D.nfrequencies, 1:Dnew.nsamples, i) = d;
    
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

D = Dnew;
D = D.history('spm_eeg_weight_epochs_TF', S);

save(D);

spm('Pointer', 'Arrow');
