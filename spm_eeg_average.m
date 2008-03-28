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
% $Id: spm_eeg_average.m 1278 2008-03-28 18:38:11Z stefan $

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
    error(sprintf('Trouble reading file %s', D));
end

% try
%     circularise = S.circularise_phase;
% catch
%         circularise = 0;
%     try
%        phs = D.phs;
%            if phs
%           circularise = spm_input('Straight or vector (eg PLV) average of phase angles?','+1', 'straight|vector',[0 1], 1);
%            end
%     end
% end


spm('Pointer', 'Watch'); drawnow;

% if isfield(D, 'Nfrequencies');
%     D.scale = zeros(D.Nchannels, 1, 1, D.events.Ntypes);
%
%     for i = 1:D.events.Ntypes
%       clear d
%       w = find((D.events.code == D.events.types(i)) & ~D.events.reject);
%       if ~circularise               % straight average
%         d = mean(D.data(:,:,:,w), 4);
%       else              % vector average (eg PLV for phase)
%         for c=1:D.Nchannels
%            tmp = D.data(c,:,:,w);
%            tmp = cos(tmp) + sqrt(-1)*sin(tmp);
%            d(c,:,:) = squeeze(abs(mean(tmp,4)) ./ mean(abs(tmp),4));
%         end
%       end
%       ni(i) = length(find(w));
%       D.scale(:, 1, 1, i) = max(max(abs(d), [], 3), [], 2)./32767;
%       d = int16(d./repmat(D.scale(:, 1, 1, i), [1, D.Nfrequencies, D.Nsamples]));
%       fwrite(fpd, d, 'int16');
%     end
%     Ntypes = D.events.Ntypes;
% else

% generate new meeg object with new filenames
Dnew = clone(D, ['m' fnamedat(D)], [D.nchannels D.nsamples D.nconditions]);
        
try artefact = D.artefact; catch artefact = []; end
if ~isempty(artefact)
    weights = artefact.weights;
    try thresholded = D.thresholded; catch thresholded = []; end
end
cl = unique(conditions(D));

if isfield(D, 'weights');
    d = zeros(D.nchannels, D.nsamples);

    for i = 1:D.nconditions

        for j = 1:D.nchannels
            tempwf=[];
            ti=0;
            ts=0;
            while ts==0
                ti=ti+1;
                ts=(j==thresholded{ti});

            end
            w = intersect(pickconditions(D, deblank(cl(i,:))), find(~D.reject))';
            ni(i) = length(w);
            if isempty(ts)
                ind = pickconditions(D, deblank(cl(i,:)));
                data = squeeze(D(j, :, ind))';
                for nl = ind
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
        Dnew = putdata(Dnew, j, 1:Dnew.nsamples, ind, d);

    end
else


    spm_progress_bar('Init', D.nconditions, 'Averages done'); drawnow;
    if D.nconditions > 100, Ibar = floor(linspace(1, D.nconditions, 100));
    else Ibar = [1:D.nconditions]; end

    for i = 1:D.nconditions

        w = intersect(pickconditions(D, deblank(cl(i,:))), find(~D.reject))';
        c = zeros(1, D.ntrials);
        c(w) = 1;
        
        ni(i) = length(w);

        d = zeros(D.nchannels, D.nsamples);

        if ni(i) == 0
            warning('%s: No trials for trial type %d', D.fname, conditionlabels(D, i));
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

cl = unique(D.conditions);

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

save(Dnew);

spm('Pointer', 'Arrow');
