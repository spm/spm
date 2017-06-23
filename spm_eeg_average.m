function D = spm_eeg_average(S)
% Average each channel over trials or trial types
% FORMAT D = spm_eeg_average(S)
%
% S        - optional input struct
%    fields of S:
% D        - MEEG object or filename of M/EEG mat-file with epoched data
% S.robust      - (optional) - use robust averaging
%                 .savew  - save the weights in an additional dataset
%                 .bycondition - compute the weights by condition (1,
%                                default) or from all trials (0)
%                 .ks     - offset of the weighting function (default: 3)
% S.prefix     - prefix for the output file (default - 'm')
%
% Output:
% D        - MEEG object (also written on disk)
%__________________________________________________________________________
% Copyright (C) 2008-2017 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_average.m 7125 2017-06-23 09:49:29Z guillaume $

SVNrev = '$Rev: 7125 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG averaging'); spm('Pointer','Watch');

%-Get input parameters
%--------------------------------------------------------------------------
if ~isfield(S, 'prefix'),       S.prefix = 'm';           end
if ~isfield(S, 'robust'),       S.robust = 0;             end

%-Configure robust averaging
%--------------------------------------------------------------------------
if isstruct(S.robust)
    if ~isfield(S.robust, 'savew'),        S.robust.savew =  0;        end
    if ~isfield(S.robust, 'bycondition'),  S.robust.bycondition = 0;   end
    if ~isfield(S.robust, 'ks'),           S.robust.ks =  3;           end
    if ~isfield(S.robust, 'removebad'),    S.robust.removebad =  0;    end
    
    robust      = 1;
    savew       = S.robust.savew;
    bycondition = S.robust.bycondition;
    ks          = S.robust.ks;
    removebad   = S.robust.removebad;
else
    robust = 0;
end

%-Get MEEG object
%--------------------------------------------------------------------------
D = spm_eeg_load(S.D);

montage_ind = D.montage('getindex');
D           = D.montage('switch', 0);

%-Check that there is any good data available
%--------------------------------------------------------------------------
if ntrials(D)==0 || isempty(indtrial(D, D.condlist, 'GOOD'))
    warning('No good trials for averaging were found. Nothing to do.');
    return;
end

%-Redirect to Time-Frequency averaging if necessary
%--------------------------------------------------------------------------
if strncmpi(D.transformtype,'TF',2) % TF and TFphase
    D = spm_eeg_average_TF(S);
    return
end

%-Generate new MEEG object with new files
%--------------------------------------------------------------------------
Dnew = clone(D, [S.prefix fname(D)], [D.nchannels D.nsamples D.nconditions]);
Dnew = type(Dnew, 'evoked');

if robust && savew
    Dw = clone(D, ['W' fname(D)]);
end

%-Do the averaging
%--------------------------------------------------------------------------
cl   = D.condlist;

ni = zeros(1, D.nconditions);
for i = 1:D.nconditions
    w = indtrial(D, deblank(cl{i}), 'GOOD');
    ni(i) = length(w);
    if ni(i) == 0
        warning('%s: No trials for trial type %d', D.fname, cl{i});
    end
end

goodtrials  = indtrial(D, cl, 'GOOD');
chanind     = D.indchantype({'MEEG', 'LFP'});

if prod(size(D))*8 < spm('memory')

    %-Average everything at once if there is enough memory
    %======================================================================
    spm_progress_bar('Init', D.nconditions, 'Averages done');
    if D.nconditions > 100, Ibar = floor(linspace(1, D.nconditions, 100));
    else Ibar = 1:D.nconditions; end

    if robust && ~bycondition
        W1      = ones(D.nchannels, D.nsamples, length(goodtrials));
        Y       = D(chanind, :, goodtrials);
        if removebad
            bad    = badsamples(D, chanind, ':', goodtrials);
            Y(bad) = NaN;
        end
        [Y, W2] = spm_robust_average(Y, 3, ks);
        W1(chanind, :, :) = W2;
        if savew
            Dw(:, :, goodtrials) = W1;
        end
        W = zeros(size(D));
        W(:, :, goodtrials) = W1;
    end

    for i = 1:D.nconditions

        w = indtrial(D, deblank(cl{i}), 'GOOD');

        if isempty(w)
            continue;
        end

        if ~robust
            Dnew(:, :, i) = mean(D(:, :, w), 3);
        else
            if bycondition
                W       = ones(D.nchannels, D.nsamples, length(w));
                Y       = D(chanind, :, w);
                if removebad
                    bad     = badsamples(D, chanind, ':', w);
                    Y(bad)  = NaN;
                end
                [Y, W1] = spm_robust_average(Y, 3, ks);
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
                Y       = D(j, :, goodtrials);
                if removebad
                    bad     = badsamples(D, j, ':', goodtrials);
                    Y(bad)  = NaN;
                end
                [Y, W1] = spm_robust_average(Y, 3, ks);
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

            w = indtrial(D, deblank(cl{i}), 'GOOD');

            if isempty(w)
                continue;
            end

            if ~robust || ~ismember(j, chanind)
                Dnew(j, :, i) = mean(D(j, :, w), 3);
            else
                if bycondition
                    Y       = D(j, :, w);
                    if removebad
                        bad     = badsamples(D, j, ':', w);
                        Y(bad)  = NaN;
                    end
                    [Y, W] = spm_robust_average(Y, 3, ks);
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
Dnew = conditions(Dnew, ':', cl);
Dnew = repl(Dnew, ':', ni);
Dnew = montage(Dnew, 'switch', montage_ind);

%-Display averaging statistics
%--------------------------------------------------------------------------
fprintf('%s: Number of replications per contrast:\n', Dnew.fname);      %-#
for i = 1:D.nconditions
    fprintf('  average %s: %d trials\n', cl{i}, ni(i));                 %-#
end

if robust
    disp('Robust averaging might have introduced high frequencies in the data. It is advised to re-apply low-pass filter.');
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

%-Cleanup
%--------------------------------------------------------------------------
fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
spm('FigName','M/EEG averaging: done'); spm('Pointer','Arrow');
