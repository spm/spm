function D = spm_eeg_rereference(S);
% rereference EEG data to new reference channel(s)
% FORMAT D = spm_eeg_rereference(S)
%
% S		    - struct (optional)
% (optional) fields of S:
% D			- filename of EEG-file with continuous data
% Nnewref   - number of new reference channels
% newref    - new reference channels
%
% Output:
% D			- EEG data struct (also written to files)
%_______________________________________________________________________
% 
% spm_eeg_rereference works under several assumptions. These are:
% (i) We assume only one channel as old reference channel
% (ii) We assume that the user never wants to 're-re-reference', but rather
%      would re-do the first rereferencing.
% (iii) The single reference be replaced by single or multiple references.
% (iv) If the new reference is only one channel, that channel is deleted from the
%       data. Multiple channels are kept in the data.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id$


[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG rereference setup',0);

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

D.fnamedat = ['R' D.fnamedat];

fpd = fopen(fullfile(P, D.fnamedat), 'w');

if ~isfield(D, 'reref')
	D.reref = [];
end

if ~isfield(D.reref, 'Nnewref')
	D.reref.Nnewref = ...
	 	spm_input('How many new reference channels', '+1', 'i', '', 1);
end

if ~isfield(D.reref, 'newref')
	for i = 1 : D.reref.Nnewref
		D.reref.newref{i} = ...
			spm_input(sprintf('New reference channel %d', i), '+1', 's', '');
	end
end

% Identify new reference channel(s)
for i = 1:D.reref.Nnewref
	tmp = find(strcmpi(D.reref.newref{i}, D.channels.name));
	if isempty(tmp)
		error('Could not find new reference %s in data set',...
							D.reref.newref{i});
	elseif tmp == D.channels.heog | tmp == D.channels.veog
		error('Cannot include EOG as reference');
    else
        Cnewref(i) = tmp;
	end
end

% get time course of new reference for all events
Tref = zeros(D.Nsamples, D.Nevents);
for i = 1:D.Nevents
	if D.reref.Nnewref > 1
		Tref(:,i) = mean(squeeze(D.data(Cnewref, :, i)))';
	else
		Tref(:,i) = squeeze(D.data(Cnewref, :, i))';
	end
end

Csetup = load(fullfile(spm('dir'), 'templates', D.channels.ctf));

in_no_eog = setdiff([1:D.Nchannels], [D.channels.heog D.channels.veog]);

oldref = [];
if isfield(D.channels, 'reference') & ~isempty(D.channels.reference)
    try
        oldref = D.channels.order_ref;
    catch 
        oldref = find(strcmpi(D.channels.reference, Csetup.Cnames));
    end
end

if isempty(oldref)
    warning(sprintf('Could not find old reference in channel template file'));
else        
    % add old reference channel
    D.channels.order(D.Nchannels+1) = oldref;
    D.channels.name{D.Nchannels+1} = D.channels.reference;
    
    % take the mean of all other (non-EOG) channels for calibration etc. 
    D.channels.base(D.Nchannels+1) = mean(D.channels.base(in_no_eog));
    D.channels.sens(D.Nchannels+1) = mean(D.channels.sens(in_no_eog));
    D.channels.calib(D.Nchannels+1) = mean(D.channels.calib(in_no_eog));
    D.Nchannels = D.Nchannels + 1;
end

% save the data
fpd = fopen(fullfile(P, D.fnamedat), 'w');

if D.reref.Nnewref == 1
    D.Nchannels = D.Nchannels - 1;
end
d = zeros(D.Nchannels, D.Nsamples);
D.scale.dim = [1 3];
D.scale.values = zeros(D.Nchannels, D.Nevents);

% Check in new reference channels
D.channels.reference = {};
D.channels.order_ref = [];
for i = 1:length(Cnewref)
    D.channels.reference{i} = D.channels.name{Cnewref(i)};
    D.channels.order_ref(i) = Cnewref(i);
end

% delete new reference channel from channel struct
if D.reref.Nnewref == 1
    D.channels.order(Cnewref) = [];
    D.channels.name(Cnewref) = [];
    D.channels.base(Cnewref) = [];
    D.channels.sens(Cnewref) = [];
    D.channels.calib(Cnewref) = [];
end

    
for i = 1:D.Nevents
    
    % concatenate data for each event
    d = squeeze(D.data(:, :, i));
    
    % Subtract new reference time series from other channels
    d(in_no_eog,:) = d(in_no_eog,:) - repmat(Tref(:,i)', length(in_no_eog), 1);
    
    % remove data of new (single) reference channel
    if D.reref.Nnewref == 1
        d(Cnewref, :) = [];
    end

    if ~isempty(oldref)
        % data for the old reference channel
        d(end+1, :) = - Tref(:,i)';
    end
    
    
	D.scale.values(:, i) = max(abs(d'))./32767;
	if ~all(D.scale.values(:, i) == 0)
		d = int16(d./repmat(D.scale.values(:, i), 1, D.Nsamples));
	end
	
	fwrite(fpd, d, 'int16');
						
end
fclose(fpd);

D.data = [];

D.fname = ['R' D.fname];
D.datatype = 'int16';

if str2num(version('-release'))>=14
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

spm('Pointer', 'Arrow');


