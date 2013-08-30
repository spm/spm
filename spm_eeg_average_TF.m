function D = spm_eeg_average_TF(S)
% Average each channel over trials or trial types, for time-frequency data
% FORMAT D = spm_eeg_average_TF(S)
%
% S         - optional input struct
%    fields of S:
% S.D           - MEEG object or filename of M/EEG mat-file with epoched TF data
% S.circularise - flag that indicates whether average is straight (0) or
%                 vector (1) of phase angles.
% S.robust      - (optional) - use robust averaging (only for power)
%                 .savew  - save the weights in an additional dataset
%                 .bycondition - compute the weights by condition (1,
%                                default) or from all trials (0)
%                 .ks     - offset of the weighting function (default: 3)
%
% Output:
% D         - MEEG object (also written to disk).
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_average_TF.m 5624 2013-08-30 11:06:38Z vladimir $

SVNrev = '$Rev: 5624 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG TF averaging'); spm('Pointer','Watch');


if ~isfield(S, 'circularise'),       S.circularise = 0;          end

%-Configure robust averaging
%--------------------------------------------------------------------------
if isstruct(S.robust)
    if ~isfield(S.robust, 'savew'),        S.robust.savew =  0;        end
    if ~isfield(S.robust, 'bycondition'),  S.robust.bycondition = 0;   end
    if ~isfield(S.robust, 'ks'),           S.robust.ks =  3;           end
    
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

%-Check data type
%--------------------------------------------------------------------------
if ~strcmp(D.type, 'single')
    error('This function can only be applied to single trial data');
end

if ~strncmp(D.transformtype, 'TF', 2) % TF and TFphase
    error('This function can only be applied to time-frequency data.');
end

if strcmp(D.transformtype, 'TFphase')
    if ~isequal(S.robust, 0)
        warning('Robust averaging is not applicable to phase data and will not be used.');
        S.robust = 0;
    end
    plv = S.circularise;
end

%-Generate new MEEG object with new files
%--------------------------------------------------------------------------
Dnew = clone(D, [S.prefix fname(D)], [D.nchannels D.nfrequencies D.nsamples D.nconditions]);

if robust && savew
    Dw = clone(D, ['W' fname(D)]);
end

cl   = D.condlist;

spm_progress_bar('Init', D.nchannels, 'Channels completed');

ni = zeros(1,D.nconditions);
for i = 1:D.nconditions
    w = indtrial(D, deblank(cl{i}), 'GOOD');
    ni(i) = length(w);
    if ni(i) == 0
        warning('%s: No trials for trial type %s', D.fname, cl{i});
    end
end

goodtrials  =  indtrial(D, cl, 'GOOD');

for j = 1:D.nchannels
    if robust && ~bycondition
        Y       = D(j, :, :, goodtrials);
        if removebad
            bad     = badsamples(D, j, ':', goodtrials);
            bad     = reshape(bad, [size(bad, 1), 1, size(bad, 2), size(bad, 3)]);
            bad     = repmat(bad, [1 D.nfrequencies, 1, 1]);
            Y(bad)  = NaN;
        end
        [Y, W1] = spm_robust_average(Y, 4, ks);
        if savew
            Dw(j, :, :, goodtrials) = W1;
        end
        W = zeros([1 D.nfrequencies D.nsamples D.ntrials]);
        W(1, :, :, goodtrials) = W1;
    end
    for i = 1:D.nconditions
        
        w = indtrial(D, deblank(cl{i}), 'GOOD');
        
        if isempty(w)
            continue;
        end
        
        %-Straight average
        %------------------------------------------------------------------
        if ~strcmp(D.transformtype, 'TFphase')
            if ~robust
                Dnew(j, :, :, i) = mean(D(j, :, :, w), 4);
            else
                if bycondition
                    Y      = D(j, :, :, w);
                    if removebad
                        bad     = badsamples(D, j, ':', w);
                        bad     = reshape(bad, [size(bad, 1), 1, size(bad, 2), size(bad, 3)]);
                        bad     = repmat(bad, [1 D.nfrequencies, 1, 1]);
                        Y(bad)  = NaN;
                    end
                    [Y, W] = spm_robust_average(Y, 4, ks);
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
            tmp = exp(sqrt(-1)*tmp);
            if plv
                Dnew(j, :, :, i) = abs(mean(tmp,4));
            else
                Dnew(j, :, :, i) = angle(mean(tmp,4));
            end
        end
    end
    spm_progress_bar('Set', j);
end

spm_progress_bar('Clear');

Dnew = type(Dnew, 'evoked');

%-Update some header information
%--------------------------------------------------------------------------
Dnew = conditions(Dnew, ':', cl);
Dnew = repl(Dnew, ':', ni);

%-Display averaging statistics
%--------------------------------------------------------------------------
fprintf('%s: Number of replications per contrast:', Dnew.fname);  %-#

s = [];
for i = 1:D.nconditions
    s = [s sprintf('average %s: %d trials', cl{i}, ni(i))];
    if i < D.nconditions
        s = [s ', '];
    else
        s = [s '\n'];
    end
end
fprintf(s);                                                       %-#

%-Save new evoked M/EEG dataset
%--------------------------------------------------------------------------
Dnew = Dnew.history('spm_eeg_average', S);
save(Dnew);

if robust && savew
    Dw = Dw.history('spm_eeg_average', S);
    save(Dw);
end

D = Dnew;

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG TF averaging: done'); spm('Pointer', 'Arrow');
