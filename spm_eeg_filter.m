function D = spm_eeg_filter(S)
% low-pass filter EEG data
% FORMAT D = spm_eeg_filter(S)
% 
% S       - struct (optional)
% (optional) fields of S:
% D		  - filename of EEG-data file or EEG data struct
% filter  - struct with the following fields:
%    type		- type of filter, currently only 'butterworth'
%    PHz        - cutoff [Hz]
%    parameter	- filter coefficients
% 
% D			- EEG data struct (also written to files)
%_______________________________________________________________________
% 
% spm_eeg_filter low-pass filters EEG/MEG epoched data.
%_______________________________________________________________________
% @(#)spm_eeg_filter.m	1.1 Stefan Kiebel 04/06/28

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG lowpass filter setup',0);

try
    D = S.D;
catch
    D = spm_get(1, '.mat', 'Select EEG mat file');
end

P = spm_str_manip(D, 'H');

try
	D = spm_eeg_ldata(D);
catch    
	error(sprintf('Trouble reading file %s', D));
end
	
if ~isfield(S.filter, 'type')
	D.filter.type =...
	 	spm_input('filter type', '+1', 'b', 'butterworth', [1], 1);
end

try
    PHz = S.filter.PHz;
catch
    
    str = 'Cutoff [Hz]';
    YPos = -1;
    while 1
        if YPos == -1
            YPos = '+1';
        end
        [PHz, YPos] = spm_input(str, YPos, 'r');
        if PHz > 0 & PHz < D.Radc, break, end
        str = 'Cutoff must be > 0 & < sample rate';
    end
    
    if D.filter.type == 1
        % butterworth
        D.filter.para = [];
        [B, A] = butter(min(5, round(D.Nsamples/20)), PHz/D.Radc);
        D.filter.para(1, :) = B;
        D.filter.para(2, :) = A;
    end
end

spm('Pointer', 'Watch');

% Prepare for writing data
D.fnamedat = ['f' D.fnamedat];
fpd = fopen(fullfile(P, D.fnamedat), 'w');

d = zeros(D.Nchannels, D.Nsamples);

% Assemble epochs
k = 1;
D.scale.dim = [1 3];
D.scale.values = zeros(D.Nchannels, D.Nevents);

for i = 1:D.Nevents

	d = squeeze(D.data(:, :, i));

	for j = 1:D.Nchannels
		d(j,:) = filter(B, A, d(j,:));
	end
    
	D.scale.values(:, k) = max(abs(d'))./32767;
	if ~all(D.scale.values(:, k) == 0)
		d = int16(d./repmat(D.scale.values(:, k), 1, D.Nsamples));
	end
	fwrite(fpd, d, 'int16');
						
	index(k) = i;
	k = k +1;
end

fclose(fpd);

D.data = [];

D.fname = ['f' D.fname];
 
save(fullfile(P, D.fname), 'D');

spm('Pointer', 'Arrow');
