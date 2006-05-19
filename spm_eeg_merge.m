function Dout = spm_eeg_merge(S);
% FORMAT D = spm_eeg_merge(S)
%
% S		    - optional input struct
% (optional) fields of S:
% D			- filename of EEG mat-file with continuous data
% recode    - a cell vector with one cell for each file. A vector in each cell
%             codes the new event type. This allows to either recode events
%             or keep their old codes.
% Output:
% D			- EEG data struct (also written to files)

% concatenates epoched single trial files
% concatenated file overwrites first file selected
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_merge.m 539 2006-05-19 17:59:30Z Darren $

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

    % checks about same number of channels, Nsamples and Radc
    if D{1}.Nchannels ~= D{i}.Nchannels
        error(sprintf('Data don''t have the same number of channels.\nThere is a difference between files %s and %s.', D{1}.fname, D{i}.fname));
    end

    if D{1}.Nsamples ~= D{i}.Nsamples
        error(sprintf('Data don''t have the same number of time points.\nThere is a difference between files %s and %s.', D{1}.fname, D{i}.fname));
    end

    if D{1}.Radc ~= D{i}.Radc
        error(sprintf('Data don''t have the same sampling rate.\nThere is a difference between files %s and %s.', D{1}.fname, D{i}.fname));
    end

	Dtmp = D{i};

	Dout.Nevents = Dout.Nevents + Dtmp.Nevents;

    % recode event types!
    try
        S.recode{i};
    catch
        S.recode{i} = spm_input(sprintf('Types: %s', spm_str_manip(Dtmp.fname, 'r')),...
            '+1', 'n', sprintf('%d ', unique(Dtmp.events.code)), Dtmp.events.Ntypes)';
    end

    if ~isempty(Dtmp.events.code)
        % recode labels
        code_new = Dtmp.events.code;
        for j = 1:Dtmp.events.Ntypes
            code_new(find(Dtmp.events.code == Dtmp.events.types(j))) = S.recode{i}(j);
        end
    else
        code_new = S.recode{i};
    end
    
    Dout.events.code = [Dout.events.code code_new];
        
    Dout.events.time = [Dout.events.time [Dtmp.events.time] + Dtmp.Nsamples];
        
	if isfield(Dout.events, 'reject')
		Dout.events.reject = [Dout.events.reject Dtmp.events.reject];
    end
    
    Dout.channels.Bad = unique([Dout.channels.Bad Dtmp.channels.Bad]);
end

Dout.events.types = unique(Dout.events.code);
Dout.events.Ntypes = length(Dout.events.types);

fpd = fopen(fullfile(P, Dout.fnamedat), 'w');

Dout.scale = zeros(Dout.Nchannels, 1, Dout.Nevents);

k = 1;
for i = 1:Nfiles
    for j = 1:D{i}.Nevents

        d = squeeze(D{i}.data(:, :, j));
        Dout.scale(:, 1, k) = spm_eeg_write(fpd, d, 2, Dout.datatype);
        k = k + 1;
    end
end

fclose(fpd);

Dout.data = [];

D = Dout;

if spm_matlab_version_chk('7') >= 0
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

spm('Pointer', 'Arrow');



