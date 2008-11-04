function D = spm_eeg_average_TF(S)
% averages each channel over trials or trial types, for time-frequency data.
% FORMAT D = spm_eeg_average_TF(S)
%
% S         - optional input struct
% (optional) fields of S:
% D                 - filename of EEG mat-file with epoched data
% circularise       - flat that indicates whether average is straight (0) or
%                     vector (1) of phase angles.
% Output:
% D         - EEG data struct (also written to files)
%_______________________________________________________________________
%
% spm_eeg_average_TF averages single trial time-frequency data within trial type.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_average_TF.m 2438 2008-11-04 11:21:19Z stefan $

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
    error(sprintf('Trouble reading file %s', D));
end

try
    circularise = S.circularise;
catch
    circularise = 0;
    try
        
        if strcmp(D.transformtype, 'TFphase')
            circularise = spm_input('Straight or vector (eg PLV) average of phase angles?','+1', 'straight|vector',[0 1], 1);
        end
    end
    S.circularise = circularise;
end

spm('Pointer', 'Watch'); drawnow;

if ~strcmp(D.type, 'single')
    error('This function can only be applied to single trial data');
else
    if strncmp(D.transformtype, 'TF', 2) % TF and TFphase

        % generate new meeg object with new filenames
        Dnew = clone(D, ['m' fnamedat(D)], [D.nchannels D.nfrequencies D.nsamples D.nconditions]);
        cl = unique(conditions(D));

        spm_progress_bar('Init', D.nconditions, 'Averages done'); drawnow;
        if D.nconditions > 100, Ibar = floor(linspace(1, D.nconditions, 100));
        else Ibar = [1:D.nconditions]; end

        for i = 1:D.nconditions

            w = pickconditions(D, deblank(cl{i}), 1)';

            ni(i) = length(w);

            if ni(i) == 0
                warning('%s: No trials for trial type %d', D.fname, conditionlabels(D, i));
            else
                if ~circularise % straight average

                    for j = 1:D.nchannels
                        Dnew(j, 1:Dnew.nfrequencies, 1:Dnew.nsamples, i) = mean(D(j, :, :, w), 4);
                    end
                else% vector average (eg PLV for phase)

                    for j = 1:D.nchannels
                        tmp = D(j,:,:,w);
                        tmp = cos(tmp) + sqrt(-1)*sin(tmp);
                        Dnew(j, 1:Dnew.nsamples, i) = squeeze(abs(mean(tmp,4)) ./ mean(abs(tmp),4));
                    end
                end

                if ismember(i, Ibar)
                    spm_progress_bar('Set', i);
                    drawnow;
                end
            end
        end
    else
        error('This function can only be applied to time-frequency data.');
    end
end
spm_progress_bar('Clear');

Dnew = type(Dnew, 'evoked');

% jump outside methods to reorganise trial structure
sD = struct(Dnew);

[sD.trials.label] = deal([]);
for i = 1:D.nconditions
    sD.trials(i).label = cl{i};
end

sD.trials = sD.trials(1:D.nconditions);

for i = 1:D.nconditions, sD.trials(i).repl = ni(i); end

Dnew = meeg(sD);

cl = D.condlist;

disp(sprintf('%s: Number of replications per contrast:', Dnew.fname))
s = [];
for i = 1:D.nconditions
    s = [s sprintf('average %s: %d trials', cl{i}, ni(i))];
    if i < D.nconditions
        s = [s sprintf(', ')];
    else
        s = [s '\n'];
    end
end
disp(sprintf(s))

D = Dnew;
D = D.history('spm_eeg_average_TF', S);

save(D);

spm('Pointer', 'Arrow');
