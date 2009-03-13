function D = spm_eeg_weight_epochs_TF(S)
% Compute contrasts over trials or trial types, for time-frequency data
% FORMAT D = spm_eeg_weight_epochs_TF(S)
%
% S             - optional input struct
% (optional) fields of S:
%   S.D         - MEEG object or filename of M/EEG mat-file with TF data
%   S.c         - contrast matrix, each row computes a contrast of the data
%   S.WeightAve - flag whether average should be weighted by number of
%                 replications (yes (1), no (0))
%   S.label     - a cell array of contrast labels (default: prompt)
%
% D             - MEEG object (also written to disk)
%__________________________________________________________________________
%
% spm_eeg_weight_epochs_TF computes contrasts of data, over epochs of data, 
% for time-frequency data. The input is a single MEEG file.
% The argument c must have dimensions Ncontrasts X Nepochs, where Ncontrasts
% is the number of contrasts and Nepochs the number of epochs, i.e. each 
% row of c contains one contrast vector. The output is a MEEG file with 
% Ncontrasts epochs. The typical use is to compute, for display purposes, 
% contrasts like the difference or interaction between trial types in 
% channel space. Another possible use is remove trials from the data file, 
% by using a contrast that contains zeros for the to be removed file.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_weight_epochs_TF.m 2874 2009-03-13 11:52:24Z guillaume $

SVNrev = '$Rev: 2874 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG weigth epochs TF'); spm('Pointer','Watch');

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

if ~strncmp(D.transformtype, 'TF', 2)
    error('Use this function for time-frequency data only');
end

%-Get parameters
%--------------------------------------------------------------------------
try
    c = S.c;
catch
    c = spm_input('Enter contrasts', '1', 'x', '', inf, eye(D.ntrials));
    S.c = c;
end

Ncontrasts = size(c, 1);
if ~isfield(S, 'label') || numel(S.label)~=size(c, 1)
    label = cell(1,Ncontrasts);
    for i = 1:Ncontrasts
        label{i} = spm_input(['Label of contrast ' num2str(i)], '+1', 's');
    end
    S.label = label;
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

Ncontrasts = size(c, 1);

% generate new meeg object with new filenames
Dnew = clone(D, ['m' fnamedat(D)], [D.nchannels D.nfrequencies D.nsamples Ncontrasts]);

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

    d = zeros(D.nchannels, D.nfrequencies, D.nsamples);

    for j = 1:D.nchannels
        for f = 1:D.nfrequencies
            d(j, f, :) = c(i,:) * squeeze(D(j, f, :, :))';
        end
    end

    Dnew(1:Dnew.nchannels, 1:D.nfrequencies, 1:Dnew.nsamples, i) = d;
    
    newrepl(i) = sum(D.repl(find(c(i,:)~=0)));

    if ismember(i, Ibar), spm_progress_bar('Set', i); end

end

spm_progress_bar('Clear');

Dnew = conditions(Dnew, [], S.label);
Dnew = trialonset(Dnew, [], []);
Dnew = reject(Dnew, [], 0);
Dnew = repl(Dnew, [], newrepl);

% remove previous source reconsructions
if isfield(Dnew,'inv')
    Dnew = rmfield(Dnew,'inv');
end

D = Dnew;
D = D.history('spm_eeg_weight_epochs_TF', S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG weight epochs TF: done'); spm('Pointer','Arrow');
