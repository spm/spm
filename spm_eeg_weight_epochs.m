function D = spm_eeg_weight_epochs(S)
% Compute contrasts over trials or trial types
% FORMAT D = spm_eeg_weight_epochs(S)
%
% S         - optional input struct
% (optional) fields of S:
% D         - filename of EEG mat-file with epoched data
% c         - contrast matrix, each row computes a contrast of the data
% label     - cell array of labels for the contrasts, the same size as
%             number of rows in c
% WeightAve - flag whether average should be weighted by number of
%             replications (yes (1), no (0))
% Output:
% D         - EEG data struct (also written to disk)
%__________________________________________________________________________
%
% spm_eeg_weight_epochs computes contrasts of data, over epochs of data. The
% input is a single MEEG file.
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

% Stefan Kiebel, Rik Henson
% $Id: spm_eeg_weight_epochs.m 3895 2010-05-20 11:43:57Z vladimir $

SVNrev = '$Rev: 3895 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG Contrasts'); spm('Pointer','Watch');

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

%-Get parameters
%--------------------------------------------------------------------------
if ~(isfield(S, 'c') && isfield(S, 'label') && numel(S.label)==size(S.c, 1))
    S.c     = [];
    S.label = {};
    i = 1;
    while 1
        S.c(i, :) = spm_input(['Enter contrast ' num2str(i)],  1, 'x', '', 1, eye(D.ntrials));
        S.label{i} = spm_input(['Label of contrast ' num2str(i)], '+1', 's');
        
        if spm_input('Add another?', '+1', 'yes|no', [0 1]);
            break;
        end
        i = i+1;
    end
end

c          = S.c;
Ncontrasts = size(c, 1);

% Pad with zeros as in the contrast manager
if size(c, 2) < D.ntrials
    c = [c zeros(Ncontrasts, D.ntrials - size(c, 2))];
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

if strncmp(D.transformtype, 'TF', 2)
    Dnew = clone(D, ['w' fnamedat(D)], [D.nchannels D.nfrequencies D.nsamples Ncontrasts]);
else
    % generate new meeg object with new filenames
    Dnew = clone(D, ['w' fnamedat(D)],[D.nchannels D.nsamples Ncontrasts]);
end

spm_progress_bar('Init', Ncontrasts, 'Contrasts computed');
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

    if strncmp(D.transformtype, 'TF', 2)
        d = zeros(D.nchannels, D.nfrequencies, D.nsamples);

        for j = 1:D.nchannels
            for f = 1:D.nfrequencies
                d(j, f, :) = c(i,:) * squeeze(D(j, f, :, :))';
            end
        end

        Dnew(1:Dnew.nchannels, 1:D.nfrequencies, 1:Dnew.nsamples, i) = d;

    else

        d = zeros(D.nchannels, D.nsamples);

        for j = 1:D.nchannels
            d(j, :) = c(i,:) * squeeze(D(j, :, :))';
        end

        Dnew(1:Dnew.nchannels, 1:Dnew.nsamples, i) = d;
    end
    
    newrepl(i) = sum(D.repl(find(c(i,:)~=0)));

    if ismember(i, Ibar), spm_progress_bar('Set', i); end

end

spm_progress_bar('Clear');

%-
%--------------------------------------------------------------------------
Dnew = conditions(Dnew, [], S.label);
Dnew = trialonset(Dnew, [], []);
Dnew = reject(Dnew, [], 0);
Dnew = repl(Dnew, [], newrepl);

% remove previous source reconsructions
if isfield(Dnew,'inv')
    Dnew = rmfield(Dnew,'inv');
end

%-Save new M/EEG data
%--------------------------------------------------------------------------
D = Dnew;
D = D.history('spm_eeg_weight_epochs', S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG Contrasts: done'); spm('Pointer','Arrow');
