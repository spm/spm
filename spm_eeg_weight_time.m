function D = spm_eeg_weight_time(S)
% function used for pre-multiplying epoched data with a user-specified
% time-specific weighting function
% FORMAT D = spm_eeg_weight_time(S)
% 
% S		    - optional input struct
% (optional) fields of S:
% D			- filename of EEG mat-file
% weight    - weigthing function
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_weight_time.m 539 2006-05-19 17:59:30Z Darren $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG weighting setup',0);

try
    D = S.D;
catch
    D = spm_select(1, 'mat', 'Select EEG mat file');
end

P = spm_str_manip(D, 'H');

try
	D = spm_eeg_ldata(D);
catch    
	error(sprintf('Trouble reading file %s', D));
end

try
    weight = S.weight;
catch
        weight = spm_input('Weighting function', -1, 'r', '', [1 D.Nsamples]);
end

spm('Pointer', 'Watch');drawnow;

% Prepare for writing data
D.fnamedat = ['w' D.fnamedat];
fpd = fopen(fullfile(P, D.fnamedat), 'w');

D.scale = zeros(D.Nchannels, 1, D.Nevents);

spm_progress_bar('Init', D.Nevents, 'Events weighted'); drawnow;
if D.Nevents > 100, Ibar = floor(linspace(1, D.Nevents,100));
else, Ibar = [1:D.Nevents]; end

WM = repmat(weight, D.Nchannels, 1);

for i = 1:D.Nevents
    d = squeeze(D.data(:, :, i));
    d = WM.*d;
    D.scale(:, 1, i) = spm_eeg_write(fpd, d, 2, D.datatype);
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end

end

spm_progress_bar('Clear');

fclose(fpd);

D.Wtime = weight;
D.data = [];
D.fname = ['w' D.fname];

if spm_matlab_version_chk('7') >= 0,
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

spm('Pointer', 'Arrow');
