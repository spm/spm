function D = spm_eeg_average_TF(S)
% Average each channel over trials or trial types, for time-frequency data
% FORMAT D = spm_eeg_average_TF(S)
%
% S         - optional input struct
% (optional) fields of S:
% S.D           - MEEG object or filename of M/EEG mat-file with epoched TF data
% S.circularise - flag that indicates whether average is straight (0) or
%                 vector (1) of phase angles.
% Output:
% D         - MEEG object (also written to disk)
%__________________________________________________________________________
%
% This function averages single trial time-frequency data within trial type.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_average_TF.m 2850 2009-03-10 21:54:38Z guillaume $

SVNrev = '$Rev: 2850 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG TF averaging'); spm('Pointer','Watch');

%-Get MEEG object
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D = spm_eeg_load(D);

%-Check data type
%--------------------------------------------------------------------------
if ~strcmp(D.type, 'single')
    error('This function can only be applied to single trial data');
elseif ~strncmp(D.transformtype, 'TF', 2) % TF and TFphase
    error('This function can only be applied to time-frequency data.');
end

%-Get other parameters
%--------------------------------------------------------------------------
try
    circularise = S.circularise;
catch
    circularise = 0;
    if strcmp(D.transformtype, 'TFphase')
        circularise = spm_input('Average of phase angles?','+1', ...
            'straight|vector(PLV)',[0 1], 1);
    end
    S.circularise = circularise;
end

%-Generate new MEEG object with new files
%--------------------------------------------------------------------------
Dnew = clone(D, ['m' fnamedat(D)], [D.nchannels D.nfrequencies D.nsamples D.nconditions]);
cl   = D.condlist;

spm_progress_bar('Init', D.nconditions, 'Averages done');

ni = zeros(1,D.nconditions);

for i = 1:D.nconditions

    w = pickconditions(D, deblank(cl{i}), 1)';

    ni(i) = length(w);

    if ni(i) == 0
        warning('%s: No trials for trial type %d', D.fname, conditionlabels(D, i));
    else
        
        %-Straight average
        %------------------------------------------------------------------
        if ~circularise

            for j = 1:D.nchannels
                Dnew(j, 1:Dnew.nfrequencies, 1:Dnew.nsamples, i) = mean(D(j, :, :, w), 4);
            end
            
        %-Vector average (eg PLV for phase)
        %------------------------------------------------------------------
        else

            for j = 1:D.nchannels
                tmp = D(j, :, :, w);
                tmp = cos(tmp) + sqrt(-1)*sin(tmp);
                Dnew(j, 1:Dnew.nsamples, i) = squeeze(abs(mean(tmp,4)) ./ mean(abs(tmp),4));
            end
            
        end

        spm_progress_bar('Set', i);
    end
end

spm_progress_bar('Clear');

Dnew = type(Dnew, 'evoked');

%-Reorganise trial structure
%--------------------------------------------------------------------------
sD = struct(Dnew);

[sD.trials.label] = deal([]);
for i = 1:D.nconditions
    sD.trials(i).label = cl{i};
end

sD.trials = sD.trials(1:D.nconditions);

for i = 1:D.nconditions, sD.trials(i).repl = ni(i); end

Dnew = meeg(sD);

%-Display averaging statistics
%--------------------------------------------------------------------------
cl = D.condlist;
disp(sprintf('%s: Number of replications per contrast:', Dnew.fname));  %-#
s = [];
for i = 1:D.nconditions
    s = [s sprintf('average %s: %d trials', cl{i}, ni(i))];
    if i < D.nconditions
        s = [s sprintf(', ')];
    else
        s = [s '\n'];
    end
end
disp(sprintf(s));                                                       %-#

%-Save new evoked M/EEG dataset
%--------------------------------------------------------------------------
D = Dnew;
D = D.history(mfilename, S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG TF averaging: done'); spm('Pointer', 'Arrow');
