function D = spm_eeg_average(S)
% Average each channel over trials or trial types
% FORMAT D = spm_eeg_average(S)
%
% S        - optional input struct
% (optional) fields of S:
% D        - MEEG object or filename of M/EEG mat-file with epoched data
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
% $Id: spm_eeg_average.m 2715 2009-02-09 17:10:41Z vladimir $

SVNrev = '$Rev: 2715 $';

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

%-Redirect to Time-Frequency averaging if necessary
%--------------------------------------------------------------------------
try
    if D.nfrequencies > 1
        D = spm_eeg_average_TF(S);
        return
    end
end

%-Generate new MEEG object with new files
%--------------------------------------------------------------------------
Dnew = clone(D, ['m' fnamedat(D)], [D.nchannels D.nsamples D.nconditions]);
Dnew = type(Dnew, 'evoked');

%-Do the averaging
%--------------------------------------------------------------------------
cl   = D.condlist;
try artefact = D.artefact; catch artefact = []; end
if isfield(artefact, 'weights')
    
    %-Weighted averaging
    %======================================================================
    weights = artefact.weights;
    try thresholded = D.thresholded; catch thresholded = []; end
    
    d = zeros(D.nchannels, D.nsamples);

    for i = 1:D.nconditions

        for j = 1:D.nchannels
            w      = pickconditions(D, cl{i}, true)';
            ni(i)  = length(w);
            
            ti     = 0;
            ts     = 0;
            while ts==0
                ti = ti+1;
                ts = (j == thresholded{ti});
            end
            if isempty(ts)
                ind    = pickconditions(D, cl{i}, 1);
                data   = squeeze(D(j, :, ind))';
                tempwf = [];
                for nl = ind'
                    tempwf = [tempwf, weights(j, (nl-1)*D.nsamples+1:nl*D.nsamples)];
                end
                data   = data';
                tempwf = reshape(tempwf, D.nsamples, length(ind));

                for t = 1:size(data,1)
                    B(t) = sum(tempwf(t,:).*data(t,:))/sum(tempwf(t,:));
                end
                
                if isfield(artefact, 'Smoothing')
                    sm  = gausswin(artefact.Smoothing);
                    sm  = sm/sum(sm);
                    mB  = mean(B);
                    B   = conv(sm,B-mean(B));
                    B   = B(floor(artefact.Smoothing/2):end-ceil(artefact.Smoothing/2));
                    B   = B + mB;
                end
                d(j, :) = B;
            else
                d(j,:)  = zeros(1,D.nsamples);
            end
        end % for j = 1:D.nchannels
        
        Dnew(1:Dnew.nchannels, 1:Dnew.nsamples, i) = d;

    end % for i = 1:D.nconditions
    
else
    %-Averaging
    %======================================================================
    spm_progress_bar('Init', D.nconditions, 'Averages done');
    if D.nconditions > 100, Ibar = floor(linspace(1, D.nconditions, 100));
    else Ibar = [1:D.nconditions]; end

    ni = zeros(1,D.nconditions);
    for i = 1:D.nconditions

        d     = zeros(D.nchannels, D.nsamples);
        
        w     = pickconditions(D, deblank(cl{i}), true)';
        c     = zeros(1, D.ntrials);
        c(w)  = 1;

        ni(i) = length(w);
        if ni(i) == 0
            warning('%s: No trials for trial type %d', D.fname, cl{i});
        else
            c = c ./ sum(c); % vector of trial-wise weights
            for j = 1:D.nchannels
                d(j, :) = c * squeeze(D(j, :, :))';
            end
        end
        
        Dnew(1:Dnew.nchannels, 1:Dnew.nsamples, i) = d;

        if ismember(i, Ibar), spm_progress_bar('Set', i); end

    end % for i = 1:D.nconditions
end

spm_progress_bar('Clear');

%-Reorganise trial structure
%--------------------------------------------------------------------------
sD = struct(Dnew);

[sD.trials.label] = deal([]);
for i = 1:D.nconditions
    sD.trials(i).label = cl{i};
end

sD.trials = sD.trials(1:D.nconditions);

for i = 1:D.nconditions, sD.trials(i).repl = ni(i); end
if isfield(sD.other, 'artefact');
    sD.other = rmfield(sD.other, 'artefact');
end

Dnew = meeg(sD);

%-Display averaging statistics
%--------------------------------------------------------------------------
disp(sprintf('%s: Number of replications per contrast:', Dnew.fname));  %-#
s  = [];
cl = D.condlist;
for i = 1:D.nconditions
    s = [s sprintf('average %s: %d trials', cl{i}, ni(i))];
    if i < D.nconditions
        s = [s ', '];
    else
        s = [s '\n'];
    end
end
disp(sprintf(s));                                                       %-#

%-Save new evoked M/EEG dataset
%--------------------------------------------------------------------------
D = Dnew;
D = D.history('spm_eeg_average', S);
save(D);

%-Eventually display it
%--------------------------------------------------------------------------
if ~isfield(S, 'review') || S.review
    spm_eeg_review(D);
end

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG averaging: done'); spm('Pointer','Arrow');
