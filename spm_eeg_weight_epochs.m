function D = spm_eeg_weight_epochs(S);
% computes contrasts over trials or trial types.
% FORMAT D = spm_eeg_weight_epochs(S)
%
% S		    - optional input struct
% (optional) fields of S:
% D			- filename of EEG mat-file with epoched data
% c         - contrast matrix, each column computes a contrast of the data
%
% Output:
% D			- EEG data struct (also written to files)
%_______________________________________________________________________
%
% spm_eeg_weight_epochs computes contrasts of data, over epochs of data. The
% input is a single MEEG file.
% The argument c must have dimensions N_contrasts x N_epochs, where N_contrasts is
% the number of contrasts and N_epochs the number of epochs, i.e. each row of c
% contains one contrast vector. The output
% is a MEEG file with N_contrasts epochs. The typical use is to compute,
% for display purposes, contrasts like the difference or interaction
% between trial types in channel space.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_weight_epochs.m 607 2006-08-31 12:29:39Z james $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG averaging setup',0);

try
    D = S.D;
catch
    D = spm_select(1, '.*\.mat$', 'Select EEG mat file');
end

P = spm_str_manip(D, 'H');

try
    D = spm_eeg_ldata(D);
catch
    error(sprintf('Trouble reading file %s', D));
end

try
    c = S.c;
catch
    c = spm_input('Enter contrasts', '1', 'x', '', inf, eye(D.Nevents));
end

% incorporate sanity check of c, later

D.fnamedat = ['m' D.fnamedat];

fpd = fopen(fullfile(P, D.fnamedat), 'w');

spm('Pointer', 'Watch'); drawnow;

N_contrasts = size(c, 1);

D.scale = zeros(D.Nchannels, 1, N_contrasts);

spm_progress_bar('Init', N_contrasts, 'Contrasts computed'); drawnow;
if N_contrasts > 100, Ibar = floor(linspace(1, N_contrasts, 100));
else, Ibar = [1:N_contrasts]; end

for i = 1:N_contrasts

    d = zeros(D.Nchannels, D.Nsamples);

    %if D.events.repl(i) == 0
      %  warning('%s: No trials for contrast %d', D.fname, i);
      % else
        for j = 1:D.Nchannels
            d(j, :) = c(i,:)*squeeze(D.data(j, :, :))';
        end
    %end

    D.scale(:, 1, i) = spm_eeg_write(fpd, d, 2, D.datatype);

    if ismember(i, Ibar)
        spm_progress_bar('Set', i);
        drawnow;
    end

end

spm_progress_bar('Clear');

fclose(fpd);

D.Nevents = N_contrasts;
D.events.code = [1:N_contrasts];

D.events.time = [];
D.events.types = D.events.code;
D.events.Ntypes = N_contrasts;

D.data = [];
D.events.reject = zeros(1, D.Nevents);
D.events.blinks = zeros(1, D.Nevents);
D.events. repl = zeros(1, D.Nevents);

D.fname = ['m' D.fname];

if spm_matlab_version_chk('7') >= 0,
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

spm('Pointer', 'Arrow');
