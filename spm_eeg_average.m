function D = spm_eeg_average(S)
% averages each channel over trials or trial types.
% FORMAT D = spm_eeg_average(S)
%
% S         - optional input struct
% (optional) fields of S:
% D         - filename of EEG mat-file with epoched data
%
% Output:
% D         - EEG data struct (also written to files)
%_______________________________________________________________________
%
% spm_eeg_average averages single trial data within trial type.
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_average.m 2398 2008-10-24 10:32:34Z vladimir $

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


spm('Pointer', 'Watch'); drawnow;

% generate new meeg object with new filenames
Dnew = clone(D, ['m' fnamedat(D)], [D.nchannels D.nsamples D.nconditions]);
        
try artefact = D.artefact; catch artefact = []; end
if ~isempty(artefact)
    weights = artefact.weights;
    try thresholded = D.thresholded; catch thresholded = []; end
end
cl = unique(conditions(D));

if isfield(artefact, 'weights');
    d = zeros(D.nchannels, D.nsamples);

    for i = 1:D.nconditions

        for j = 1:D.nchannels
            tempwf=[];
            ti=0;
            ts=0;
            while ts==0
                ti=ti+1;
                ts = (j == thresholded{ti});

            end
            w = intersect(pickconditions(D, cl{i}), find(~reject(D)))';
            ni(i) = length(w);
            if isempty(ts)
                ind = pickconditions(D, cl{i});
                data = squeeze(D(j, :, ind))';
                for nl = ind'
                    tempwf = [tempwf, weights(j, (nl-1)*D.nsamples+1:nl*D.nsamples)];
                end
                data=data';
                tempwf = reshape(tempwf, D.nsamples, length(ind));

                for t = 1:size(data,1)
                    B(t) = sum(tempwf(t,:).*data(t,:))/sum(tempwf(t,:));
                end
                
                if isfield(artefact, 'Smoothing')
                    sm=gausswin(artefact.Smoothing);
                    sm=sm/sum(sm);
                    mB=mean(B);
                    B=conv(sm,B-mean(B));
                    B=B(floor(artefact.Smoothing/2):end-ceil(artefact.Smoothing/2));
                    B=B+mB;
                end
                d(j, :) =B;
            else
                d(j,:)=zeros(1,D.nsamples);
            end
        end
        Dnew(1:Dnew.nchannels, 1:Dnew.nsamples, i) = d;

    end
else


    spm_progress_bar('Init', D.nconditions, 'Averages done'); drawnow;
    if D.nconditions > 100, Ibar = floor(linspace(1, D.nconditions, 100));
    else Ibar = [1:D.nconditions]; end

    for i = 1:D.nconditions

        w = intersect(pickconditions(D, deblank(cl{i})), find(~D.reject))';
        c = zeros(1, D.ntrials);
        c(w) = 1;
        
        ni(i) = length(w);

        d = zeros(D.nchannels, D.nsamples);

        if ni(i) == 0
            warning('%s: No trials for trial type %d', D.fname, cl{i});
        else
            c = c./sum(c); % vector of trial-wise weights
            for j = 1:D.nchannels
                d(j, :) = c*squeeze(D(j, :, :))';
            end
        end
        Dnew(1:Dnew.nchannels, 1:Dnew.nsamples, i) = d;

        if ismember(i, Ibar)
            spm_progress_bar('Set', i);
            drawnow;
        end

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
if isfield(sD.other, 'artefact');
    sD.other=rmfield(sD.other, 'artefact');
end

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
D = D.history('spm_eeg_average', S);

save(D);

if ~isfield(S, 'review') || S.review
    spm_eeg_review(D);
end

spm('Pointer', 'Arrow');
