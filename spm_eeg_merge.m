function Dout = spm_eeg_merge(S);
% FORMAT D = spm_eeg_merge(S)
% concatenates epoched single trial files
% concatenated file overwrites first file selected
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_merge.m 161 2005-05-16 14:48:27Z stefan $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG merge',0);

try
    D = S.D;
catch
    D = spm_select([2 inf], '\.mat$', 'Select EEG mat files');
end

P = spm_str_manip(D(1,:), 'H');

try
    for i = 1:size(D, 1)
       F{i} = spm_eeg_ldata(deblank(D(i, :))); 
    end
    D = F;
catch
	error('Trouble reading files');
end

Nfiles = length(D);

if Nfiles < 2
    error('Need at least two files for merging');
end

spm('Pointer', 'Watch');

Dout = D{1};
Dout.fnamedat = ['c' spm_str_manip(D{1}.fnamedat, 't')];
Dout.fname = ['c' spm_str_manip(D{1}.fname, 't')];

for i = 2:Nfiles

    % should check: same channels, Nsamples, scale, start/stop, reref, Radc, etc...

	Dtmp = D{i};

	Dout.Nevents = Dout.Nevents + Dtmp.Nevents;

    % recode event types!
    try
        S.recode{i};
    catch
        S.recode{i} = spm_input(sprintf('Recode %s', Dtmp.fname),...
            '+1', 'n', sprintf('%d ', unique(Dtmp.events.code)), max(Dtmp.events.types))';
    end

    Dout.events.code = [Dout.events.code S.recode{i}(Dtmp.events.code)];
        
    Dout.events.time = [Dout.events.time [Dtmp.events.time] + Dtmp.Nsamples];
        
	if isfield(Dout.events, 'reject')
		Dout.events.reject = [Dout.events.reject Dtmp.events.reject];
    end
    
    Dout.channels.Bad = unique([Dout.channels.Bad Dtmp.channels.Bad]);
end

Dout.events.types = unique(Dout.events.code);
Dout.events.Ntypes = length(Dout.events.types);

fpd = fopen(fullfile(P, Dout.fnamedat), 'w');

Dout.scale.dim = [1 3];
Dout.scale.values = zeros(Dout.Nchannels, Dout.Nevents);

k = 1;
for i = 1:Nfiles
    for j = 1:D{i}.Nevents

        d = squeeze(D{i}.data(:, :, j));
        Dout.scale.values(:, k) = spm_eeg_write(fpd, d, 2, Dout.datatype);
        k = k + 1;
    end
end

fclose(fpd);

Dout.data = [];

D = Dout;

if str2num(version('-release'))>=14
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

spm('Pointer', 'Arrow');



