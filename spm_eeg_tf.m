function D = spm_eeg_tf(S)
% compute instantaneous power and phase in peri-stimulus time and frequency
% FORMAT D = spm_eeg_tf(S)
% 
% D		- filename of EEG-data file or EEG data struct
% stored in struct D.events:
% fmin			- minimum frequency
% fmax			- maximum frequency
% rm_baseline	- baseline removal (1/0) yes/no
% Mfactor       - Morlet wavelet factor (can not be accessed by GUI)
% 
% D				- EEG data struct with time-frequency data (also written to files)
%_______________________________________________________________________
%
% spm_eeg_tf estimates instantaneous power and phase of data using the
% continuous Morlet wavelet transform.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_tf.m 568 2006-07-04 17:59:56Z stefan $


[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG time-frequency setup',0);

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

% Check for arguments in D
if ~isfield(D, 'tf')
	D.tf = [];
end

try
    D.tf.frequencies = S.frequencies;
catch
    D.tf.frequencies = ...
	 	spm_input('Frequencies (Hz)', '+1', 'r', '', [1, inf]);
end

try
    D.tf.rm_baseline = S.rm_baseline;
catch
    D.tf.rm_baseline = ...
		spm_input('Removal of baseline', '+1', 'y/n', [1,0], 2);
end

if D.tf.rm_baseline
    try
        D.tf.Sbaseline = S.Sbaseline;
    catch
        D.tf.Sbaseline = ...
            spm_input('Start and stop of baseline (samples)', '+1', 'i', '', 2);
    end
end

try
    D.tf.Mfactor = S.Mfactor;
catch
    D.tf.Mfactor = ...
        spm_input('Which Morlet wavelet factor?', '+1', 'r', '7', 1);
end

% NB: D.tf.channels maps directly into the data. To retrieve the position of the channel,
% use D.channels.order
try
    D.tf.channels = S.channels;
catch
    D.tf.channels = ...
	 	spm_input('Select channels', '+1', 'i', num2str([1:D.Nchannels]));
end

spm('Pointer', 'Watch'); drawnow;

M = spm_eeg_morlet(D.tf.Mfactor, 1000/D.Radc, D.tf.frequencies);

D.Nfrequencies = length(D.tf.frequencies);
fnamedat = D.fnamedat;

D.fnamedat = ['t1_' fnamedat];
D2 = D;
D2.fnamedat = ['t2_' fnamedat];

fpd = fopen(fullfile(P, D.fnamedat), 'w');
fpd2 = fopen(fullfile(P, D2.fnamedat), 'w');

D.scale = zeros(length(D.tf.channels), 1, 1, D.Nevents);
D.datatype = 'int16';

D2.scale = zeros(length(D.tf.channels), 1, 1, D.Nevents);
D2.datatype = 'int16';

for k = 1 : D.Nevents
    
	d = zeros(length(D.tf.channels), D.Nfrequencies, D.Nsamples);
	d2 = zeros(length(D.tf.channels), D.Nfrequencies, D.Nsamples);
	for j = 1 : length(D.tf.channels)
		for i = 1 : D.Nfrequencies
			tmp = conv(squeeze(D.data(D.tf.channels(j),:,k)), M{i});
            
            % time shift to remove delay
            tmp = tmp([1:D.Nsamples] + (length(M{i})-1)/2);
			
            % power
            d(j, i, :) = tmp.*conj(tmp);
            
            % phase
            d2(j, i, :) = atan2(imag(tmp), real(tmp));

		end
	end
	
	% Remove baseline over frequencies and trials
	if D.tf.rm_baseline == 1
		d = spm_eeg_bc(D, d);
    end
	D.scale(:, 1, 1, k) = (max(abs(reshape(d, [length(D.tf.channels) D.Nfrequencies*D.Nsamples])'))./32767);
	d = int16(round(d./repmat(D.scale(:, 1, 1, k), [1 D.Nfrequencies D.Nsamples])));

    D2.scale(:, 1, 1, k) = (max(abs(reshape(d2, [length(D2.tf.channels) D2.Nfrequencies*D2.Nsamples])'))./32767);
	d2 = int16(round(d2./repmat(D2.scale(:, 1, 1, k), [1 D2.Nfrequencies D2.Nsamples])));

	fwrite(fpd, d, 'int16');	
	fwrite(fpd2, d2, 'int16');	
		
end

fclose(fpd);
fclose(fpd2);

D.data = [];
D2.data = [];

D.Nchannels = length(D.tf.channels);
D2.Nchannels = length(D2.tf.channels);

D.fname = ['t1_' D.fname];
D2.fname = ['t2_' D2.fname];
D1=D;
if spm_matlab_version_chk('7') >= 0,
    save(fullfile(P, D.fname), '-V6', 'D');
    D = D2;
    save(fullfile(P, D2.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
    D = D2;
    save(fullfile(P, D2.fname), 'D');
end
D=D1;
spm('Pointer', 'Arrow');
