function D = spm_eeg_contrast(S)
% Compute contrasts over trials or trial types
% FORMAT D = spm_eeg_contrast(S)
%
% S         - optional input struct
% fields of S:
% D         - filename of EEG mat-file with epoched data
% c         - contrast matrix, each row computes a contrast of the data
% label     - cell array of labels for the contrasts, the same size as
%             number of rows in c
% weighted  - flag whether average should be weighted by number of
%             replications (yes (1), no (0))
% prefix    - prefix for the output file (default - 'w')
%
% Output:
% D         - EEG data struct (also written to disk)
%__________________________________________________________________________
%
% spm_eeg_contrast computes contrasts of data, over epochs of data. The
% input is a single MEEG file.
% The argument c must have dimensions Ncontrasts X Nepochs, where Ncontrasts is
% the number of contrasts and Nepochs the number of epochs, i.e. each row of c
% contains one contrast vector. The output
% is a M/EEG file with Ncontrasts epochs. The typical use is to compute,
% for display purposes, contrasts like the difference or interaction
% between trial types in channel space.
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel, Rik Henson
% $Id: spm_eeg_contrast.m 6526 2015-08-20 10:28:36Z vladimir $

SVNrev = '$Rev: 6526 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG Contrast'); spm('Pointer','Watch');

%-Get MEEG object
%--------------------------------------------------------------------------
if ~isfield(S, 'prefix'),       S.prefix     = 'w';          end
if ~isfield(S, 'weighted'),     S.weighted   =  0;           end

if ~(isfield(S, 'c') && isfield(S, 'label') && numel(S.label)==size(S.c, 1))
    error('Invalid contrast specification.');
end

D = spm_eeg_load(S.D);
   

%--------------------------------------------------------------------------

c          = S.c;
Ncontrasts = size(c, 1);

% Pad with zeros as in the contrast manager
if size(c, 2) <= D.ntrials
    c = [c zeros(Ncontrasts, D.ntrials - size(c, 2))];
else
    error('The number of columns in the contrast matrix exceeds the number of trials.');
end

if ~isempty(D.repl)
    weighted = S.weighted;
else
    weighted = 0;
end

if strncmp(D.transformtype, 'TF', 2)
    Dnew = clone(D, [S.prefix fname(D)], [D.nchannels D.nfrequencies D.nsamples Ncontrasts]);
else
    Dnew = clone(D,  [S.prefix fname(D)],[D.nchannels D.nsamples Ncontrasts]);
end

spm_progress_bar('Init', Ncontrasts, 'Contrasts computed');
if Ncontrasts > 100, Ibar = floor(linspace(1, Ncontrasts, 100));
else Ibar = 1:Ncontrasts; end

for i = 1:Ncontrasts

    if weighted
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
                d(j, f, :) = c(i,:) * spm_squeeze(D(j, f, :, :), [1 2 4])';
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
Dnew = conditions(Dnew, ':', S.label);
Dnew = trialonset(Dnew, ':', []);
Dnew = trialtag(Dnew, ':', []);
Dnew = badtrials(Dnew, ':', 0);
Dnew = repl(Dnew, ':', newrepl);
Dnew = type(Dnew, 'evoked');

% remove previous source reconsructions
if isfield(Dnew,'inv')
    Dnew = rmfield(Dnew,'inv');
end

%-Save new M/EEG data
%--------------------------------------------------------------------------
D = Dnew;
D = D.history('spm_eeg_contrast', S);
save(D);

%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','M/EEG Contrast: done'); spm('Pointer','Arrow');
