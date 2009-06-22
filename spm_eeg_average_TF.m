function D = spm_eeg_average_TF(S)
% Average each channel over trials or trial types, for time-frequency data
% FORMAT D = spm_eeg_average_TF(S)
%
% S         - optional input struct
% (optional) fields of S:
% S.D           - MEEG object or filename of M/EEG mat-file with epoched TF data
% S.circularise - flag that indicates whether average is straight (0) or
%                 vector (1) of phase angles.
% S.robust      - (optional) - use robust averaging
%                 .savew  - save the weights in an additional dataset
%                 .bycondition - compute the weights by condition (1,
%                                default) or from all trials (0)
%                 .ks     - offset of the weighting function (default: 3)
%
% Output:
% D         - MEEG object (also written to disk)
%__________________________________________________________________________
%
% This function averages single trial time-frequency data within trial type.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_average_TF.m 3218 2009-06-22 15:45:53Z vladimir $

SVNrev = '$Rev: 3218 $';

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

if ~isfield(S, 'robust')
    robust = 0;
else
    robust = 1;
    if ~isfield(S.robust, 'savew')
        S.robust.savew = 0;
        savew = 0;
    else
        savew = S.robust.savew;
    end

    if ~isfield(S.robust, 'bycondition')
        S.robust.bycondition = 1;
        bycondition = 1;
    else
        bycondition = S.robust.bycondition;
    end

    if ~isfield(S.robust, 'ks')
        S.robust.ks = [];
        ks          = [];
    else
        ks          = S.robust.ks;
    end
end

%-Generate new MEEG object with new files
%--------------------------------------------------------------------------
Dnew = clone(D, ['m' fnamedat(D)], [D.nchannels D.nfrequencies D.nsamples D.nconditions]);

if robust && savew
    Dw = clone(D, ['W' fnamedat(D)], [D.nchannels D.nfrequencies D.nsamples D.ntrials]);
end

cl   = D.condlist;

spm_progress_bar('Init', D.nchannels, 'Channels completed');

ni = zeros(1,D.nconditions);
for i = 1:D.nconditions
    w = pickconditions(D, deblank(cl{i}), 1)';
    ni(i) = length(w);
    if ni(i) == 0
        warning('%s: No trials for trial type %d', D.fname, cl{i});
    end
end

goodtrials  =  pickconditions(D, cl, 1);

for j = 1:D.nchannels
    if robust && ~bycondition
        [Y, W1] = spm_robust_average(D(j, :, :, goodtrials), 4, ks);
        if savew
            Dw(j, :, :, goodtrials) = W1;
        end
        W = zeros([1 D.nfrequencies D.nsamples D.ntrials]);
        W(1, :, :, goodtrials) = W1;
    end
    for i = 1:D.nconditions

        w = pickconditions(D, deblank(cl{i}), 1)';
        
        if isempty(w)
            continue;
        end

        %-Straight average
        %------------------------------------------------------------------
        if ~circularise
            if ~robust
                Dnew(j, :, :, i) = mean(D(j, :, :, w), 4);
            else
                if bycondition
                    [Y, W] = spm_robust_average(D(j, :, :, w), 4, ks);
                    Dnew(j, :, :, i) = Y;
                    if savew
                        Dw(j, :, :, w)   = W;
                    end
                else
                    X = D(j, :, :, w);
                    X(isnan(X))      = 0;
                    Dnew(j, :, :, i) = ...
                        sum(W(1, :, :, w).*X, 4)./sum(W(1, :, :, w), 4);
                end
            end

            %-Vector average (eg PLV for phase)
            %------------------------------------------------------------------
        else
            tmp = D(j, :, :, w);
            tmp = cos(tmp) + sqrt(-1)*sin(tmp);

            if ~robust
                Dnew(j, :, i) = squeeze(abs(mean(tmp,4)) ./ mean(abs(tmp),4));
            else
                if bycondition
                    [Y, W] = spm_robust_average(tmp, 4, ks);
                    aY     = sum(W(j, :, :, :).*abs(tmp), 4)./sum(W(j, :, :, :), 4);
                    if savew
                        Dw(j, :, :, w) = W;
                    end
                else
                    tmp(isnan(tmp)) = 0;
                    Y  = sum(W(1, :, :, w).*tmp)./sum(W(1, :, :, w), 4);
                    aY = sum(W(1, :, :, w).*abs(tmp))./sum(W(1, :, :, w), 4);
                end

                Dnew(j, :, :, i) = squeeze(Y./aY);

            end
        end
    end
    spm_progress_bar('Set', j);
end

spm_progress_bar('Clear');

Dnew = type(Dnew, 'evoked');

%-Update some header information
%--------------------------------------------------------------------------
Dnew = conditions(Dnew, [], cl);
Dnew = repl(Dnew, [], ni);

%-Display averaging statistics
%--------------------------------------------------------------------------
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
Dnew = Dnew.history(mfilename, S);
save(Dnew);

if robust && savew
    Dw = Dw.history(mfilename, S);
    save(Dw);
end

D = Dnew;

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG TF averaging: done'); spm('Pointer', 'Arrow');
