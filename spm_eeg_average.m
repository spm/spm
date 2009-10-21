function D = spm_eeg_average(S)
% Average each channel over trials or trial types
% FORMAT D = spm_eeg_average(S)
%
% S        - optional input struct
% (optional) fields of S:
% D        - MEEG object or filename of M/EEG mat-file with epoched data
% S.robust      - (optional) - use robust averaging
%                 .savew  - save the weights in an additional dataset
%                 .bycondition - compute the weights by condition (1,
%                                default) or from all trials (0)
%                 .ks     - offset of the weighting function (default: 3)
% review   - review data after averaging [default: true]
%
% Output:
% D        - MEEG object (also written on disk)
%__________________________________________________________________________
%
% spm_eeg_average averages single trial data within trial type.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_average.m 3497 2009-10-21 21:54:28Z vladimir $

SVNrev = '$Rev: 3497 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG averaging'); spm('Pointer','Watch');

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

%-Check that there is any good data available
%--------------------------------------------------------------------------
if ntrials(D)==0 || all(reject(D))
    warning('No good trials for averaging were found. Nothing to do.');
    return;
end

%-Redirect to Time-Frequency averaging if necessary
%--------------------------------------------------------------------------
if strncmpi(D.transformtype,'TF',2) % TF and TFphase
    D = spm_eeg_average_TF(S);
    return
end

%-Backward compatibility
%--------------------------------------------------------------------------
persistent runonce
if isfield(D, 'artefact')
    if isempty(runonce)
        warning(['Robust averaging in spm_eeg_artefact has been deprecated. ' ...
            'Use the ''robust'' option in spm_eeg_average instead.']);
        runonce = 1;
    end
    D = spm_eeg_average5(S);
    return;
end

%-Configure robust averaging
%--------------------------------------------------------------------------
if ~isfield(S, 'robust')
    robust = spm_input('Use robust averaging?','+1','yes|no',[1 0]);
    if robust
        S.robust = [];
    else
        S.robust = false;
    end
else
    robust = isa(S.robust, 'struct');
end


if robust
    if ~isfield(S.robust, 'savew')
        S.robust.savew =  spm_input('Save weights?','+1','yes|no',[1 0]);
    end
    savew = S.robust.savew;

    if ~isfield(S.robust, 'bycondition')
        S.robust.bycondition = spm_input('Compute weights by condition?','+1','yes|no',[1 0]);
    end

    bycondition = S.robust.bycondition;

    if ~isfield(S.robust, 'ks')
        S.robust.ks =  spm_input('Offset of the weighting function', '+1', 'r', '3', 1);
    end

    ks          = S.robust.ks;

end

%-Generate new MEEG object with new files
%--------------------------------------------------------------------------
Dnew = clone(D, ['m' fnamedat(D)], [D.nchannels D.nsamples D.nconditions]);
Dnew = type(Dnew, 'evoked');

if robust && savew
    Dw = clone(D, ['W' fnamedat(D)], size(D));
end

%-Do the averaging
%--------------------------------------------------------------------------
cl   = D.condlist;

ni = zeros(1, D.nconditions);
for i = 1:D.nconditions
    w = pickconditions(D, deblank(cl{i}), 1)';
    ni(i) = length(w);
    if ni(i) == 0
        warning('%s: No trials for trial type %d', D.fname, cl{i});
    end
end

goodtrials  = pickconditions(D, cl, 1);
chanind     = meegchannels(D);

if prod(size(D))*8 < spm('memory')

    %-Average everything at once if there is enough memory
    %======================================================================
    spm_progress_bar('Init', D.nconditions, 'Averages done');
    if D.nconditions > 100, Ibar = floor(linspace(1, D.nconditions, 100));
    else Ibar = [1:D.nconditions]; end

    if robust && ~bycondition
        W1      = ones(D.nchannels, D.nsamples, length(goodtrials));
        [Y, W2] = spm_robust_average(D(chanind, :, goodtrials), 3, ks);
        W1(chanind, :, :) = W2;
        if savew
            Dw(:, :, goodtrials) = W1;
        end
        W = zeros(size(D));
        W(:, :, goodtrials) = W1;
    end

    for i = 1:D.nconditions

        w = pickconditions(D, deblank(cl{i}), true)';

        if isempty(w)
            continue;
        end

        if ~robust
            Dnew(:, :, i) = mean(D(:, :, w), 3);
        else
            if bycondition
                W        = ones(D.nchannels, D.nsamples, length(w));
                [Y, W1] = spm_robust_average(D(chanind, :, w), 3, ks);
                W(chanind, :, :)    = W1;
                Dnew(chanind, :, i) = Y;
                if length(chanind)<D.nchannels
                    Dnew(setdiff(1:D.nchannels, chanind), :, i) = mean(D(setdiff(1:D.nchannels, chanind), :, w), 3);
                end
                if savew
                    Dw(:, :, w)   = W;
                end
            else
                X = D(:, :, w);
                X(isnan(X))      = 0;
                Dnew(:, :, i) = ...
                    sum(W(:, :, w).*X, 3)./sum(W(:, :, w), 3);
            end
        end

        if ismember(i, Ibar), spm_progress_bar('Set', i); end

    end  % for i = 1:D.nconditions
else
    %-Averaging by channel
    %======================================================================
    spm_progress_bar('Init', D.nchannels, 'Channels completed');
    if D.nchannels > 100, Ibar = floor(linspace(1, D.nchannels, 100));
    else Ibar = [1:D.nchannels]; end
    for j = 1:D.nchannels                
        if robust && ~bycondition 
            if ismember(j, chanind)
                [Y, W1] = spm_robust_average(D(j, :, goodtrials), 3, ks);
                W = zeros([1 D.nsamples D.ntrials]);
                W(1, :, goodtrials) = W1;
            else
                W1 = ones(1, D.nsamples, length(goodtrials));
            end

            if savew
                Dw(j, :, goodtrials) = W1;
            end
        end

        for i = 1:D.nconditions

            w = pickconditions(D, deblank(cl{i}), true)';

            if isempty(w)
                continue;
            end

            if ~robust || ~ismember(j, chanind)
                Dnew(j, :, i) = mean(D(j, :, w), 3);
            else
                if bycondition
                    [Y, W] = spm_robust_average(D(j, :, w), 3, ks);
                    Dnew(j, :, i) = Y;
                    if savew
                        Dw(j, :, w)   = W;
                    end
                else
                    X = D(j, :, w);
                    X(isnan(X))      = 0;
                    Dnew(j, :, i) = ...
                        sum(W(1, :, w).*X, 3)./sum(W(1, :, w), 3);
                end
            end

        end  % for i = 1:D.nconditions
        if ismember(j, Ibar), spm_progress_bar('Set', j); end
    end
end

spm_progress_bar('Clear');

%-Update some header information
%--------------------------------------------------------------------------
Dnew = conditions(Dnew, [], cl);
Dnew = repl(Dnew, [], ni);

%-Display averaging statistics
%--------------------------------------------------------------------------
disp(sprintf('%s: Number of replications per contrast:', Dnew.fname));  %-#
s  = [];
for i = 1:D.nconditions
    s = [s sprintf('average %s: %d trials', cl{i}, ni(i))];
    if i < D.nconditions
        s = [s ', '];
    else
        s = [s '.'];
    end
end
disp(s);                                                       %-#

if robust
    disp('Robust averaging might have introduced high frequencies in the data. It is advised to re-apply low-pass filter');
end

%-Save new evoked M/EEG dataset
%--------------------------------------------------------------------------
Dnew = Dnew.history(mfilename, S);
save(Dnew);

if robust && savew
    Dw = Dw.history(mfilename, S);
    save(Dw);
end

D = Dnew;

%-Eventually display it
%--------------------------------------------------------------------------
if ~isfield(S, 'review') || S.review
    spm_eeg_review(D);
end

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG averaging: done'); spm('Pointer','Arrow');
