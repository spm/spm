function D = spm_eeg_rdata_bdf(S)
% converts EEG data from BDF- to SPM-format
% FORMAT D = spm_eeg_rdata_bdf(S)
% 
% S       - struct (optional)
% (optional) fields of S:
% Fdata		  - filename of bdf-file
% Fchannels   - filename of channel template
% exg_name    - cell vector that code type of exogeneous channels. Allowed
%               types are 'heog', 'veog', 'reference' and 'other' 
%_______________________________________________________________________
% 
% spm_eeg_rdata_bdf reads a continuous *.bdf file, stores everything
% in struct D and saves struct D to mat-file. The data is stored separately
% in a dat-file.
% There are calls to openbdf and readbdf, which are distributed under the
% GNU-license.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_rdata_bdf.m 701 2006-11-30 12:37:39Z james $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','read BDF data setup',0);

try
    Fdata = S.Fdata;
catch
	Fdata = spm_select(1, '\.bdf$', 'Select bdf-file');
end

try
    Fchannels = S.Fchannels;
catch
    Fchannels = spm_select(1, '\.mat$', 'Select channel template file', {}, fullfile(spm('dir'), 'EEGtemplates'));
end

D.channels.ctf = spm_str_manip(Fchannels, 't');
Csetup = load(Fchannels);
F = spm_str_manip(Fdata, 'rt');
P = spm_str_manip(Fdata, 'H');

D.fnamedat = [F '.dat'];

% read header
Hdata = openbdf(Fdata);

% build header for SPM format
%--------------------------------------------------
% find out how many channels
Cind = strmatch('Active Electrode', Hdata.Head.Transducer);
D.Nchannels = length(Cind);

% Bad channels entry
D.channels.Bad = [];

% read channel names
for i = 1:D.Nchannels
    D.channels.name{i} = deblank(Hdata.Head.Label(i,:));
end

% go through 'external' channels and ask what they are
% Cexg = find(strncmpi('exg', D.channels.name, 3));

% assume that external channels are always the 8 channels after the EEG
% channels
Neeg = 2^fix(log2(D.Nchannels));
Cexg = [Neeg+1:Neeg+8];

% print names of external channels in matlab command window
disp(sprintf('External channels:'))
for i = 1:8
   disp(sprintf('%d: %s', i, D.channels.name{Neeg+i}))
end

% HEOG, VEOG and reference
try
    Cheog = S.Cheog;
catch
    Cheog = spm_input('Heog channel(s)', 'i+1', 'i', '', inf, length(Cexg));
end

try
    Cveog = S.Cveog;
catch
    Cveog = spm_input('Veog channel(s)', 'i+1', 'i', '', inf, length(Cexg));
end

try
    Creference = S.Creference;
catch
    Creference = spm_input('Reference channel(s)', 'i+1', 'i', '', inf, length(Cexg));
end



more_exg =spm_input('Do you have any other data?', '+1', 'y/n', [1,0], 2);
no_exg=0;
if more_exg==1
    no_exg=spm_input('how many extra channels?', '+1', 'i', '', 1);
    Chexg = spm_input('Enter exg (nose or others) channel(s)', 'i+1', 'i', '', inf, length(Cexg));
end    
spm('Pointer','Watch'); drawnow;

% number of EEG channels
Neeg = D.Nchannels - length(Cexg);

% save EEG channel indices
D.channels.eeg = [1:Neeg];

if Cheog == 0
    warning(sprintf('No HEOG channel.'));
    Cheog = 0;
elseif length(Cheog) > 2
    error('Only up 2 heog channels allowed');
else
    Cheog = Cheog + Neeg;
end

if Cveog == 0
    warning(sprintf('No VEOG channel.'));
    Cveog = 0;
elseif length(Cveog) > 2
    error('Only up 2 veog channels allowed');
else
    Cveog = Cveog + Neeg;
end

if Creference == 0
    warning(sprintf('No reference channel.'));
    Creference = 0;
else
    Creference = Creference + Neeg;
end

if more_exg==1
    Chexg = Chexg + Neeg;
end
        
        
% Map name of channels (EEG only) to channel order specified in channel
% template file. EOG channels will be added further below!
for i = D.channels.eeg
	index = [];
	for j = 1:Csetup.Nchannels
		if ~isempty(find(strcmpi(D.channels.name{i}, Csetup.Cnames{j})))
			index = [index j];
		end
	end
    
	if isempty(index)
		warning(sprintf('No channel named %s found in channel template file.', D.channels.name{i}));
	else
		% take only the first found channel descriptor
		D.channels.order(i) = index(1);
	end

end

% Assume that each channel has the same sampling rate
D.Radc = Hdata.Head.SampleRate(Cind(1));

% open dat-file for writing
%---------------------------------------------------
Fout = [F '.dat'];
D.fnamedat = Fout;
fpout = fopen(fullfile(P, D.fnamedat), 'w');


counter = 0;
if Cheog ~= 0
    counter = counter + 1;
    D.channels.name{Neeg + counter} = 'HEOG';
    D.channels.heog = Neeg + counter;
    index = [];
    for j = 1:Csetup.Nchannels
        if ~isempty(find(strcmpi('HEOG', Csetup.Cnames{j})))
            index = j; break;
        end
    end
    if isempty(index)
        warning(sprintf('No HEOG channel found in channel template file.'));
    else
        % take only the first found channel descriptor
        D.channels.order(Neeg+counter) = index;
    end


end

if Cveog ~= 0
    counter = counter + 1;
    D.channels.name{Neeg + counter} = 'VEOG';
    D.channels.veog = Neeg + counter;
    index = [];
    for j = 1:Csetup.Nchannels
        if ~isempty(find(strcmpi('VEOG', Csetup.Cnames{j})))
            index = j; break;
        end
    end
    if isempty(index)
        warning(sprintf('No VEOG channel found in channel template file.'));
    else
        % take only the first found channel descriptor
        D.channels.order(Neeg+counter) = index;
    end

end


D.channels.reference = 0;
D.channels.ref_name = 'NIL';

if more_exg ==1
    D.channels.exg=[];
    for kl=1:no_exg
        counter = counter + 1;
        
        D.channels.name{Neeg + counter} = ['exg',num2str(kl)];
        D.channels.exg = [ D.channels.exg ,Neeg + counter];
        index = [];
        for j = 1:Csetup.Nchannels
            if ~isempty(find(strcmpi(['exg',num2str(kl)], Csetup.Cnames{j})))
                index = j; break;
            end
        end
        if isempty(index)
            warning(sprintf('No exg channel found in channel template file.'));
        else
            % take only the first found channel descriptor
            D.channels.order(Neeg+counter) = index;
        end
        
    end
end

D.Nchannels = Neeg + counter;
if D.Nchannels < length(D.channels.name);
    D.channels.name(D.Nchannels+1:end) = [];
end

% No scaling, represent as float32
D.scale = ones(D.Nchannels, 1, 1);

D.events.code = [];
D.events.time = [];
D.Nsamples = 0;

% loop over seconds
Nblocks = Hdata.Head.NRec;

% sometimes Hdata.Head.NRec indicates more Nblocks than are actually there.
flen = dir(Fdata);
flen = flen.bytes;
Nblocks_alt = (flen-Hdata.Head.HeadLen)/Hdata.Head.SampleRate(1)/Hdata.Head.NS/3;
if Nblocks_alt < Nblocks
    warning('Nblocks too large. True Nblocks = %f', Nblocks_alt);
    Nblocks = floor(Nblocks_alt);
end

% identify channel which stores event coding
Echannel = strmatch('Status', Hdata.Head.Label);
if length(Echannel) > 1
    Echannel = Echannel(1);
end

lidiffs = 0;
lastevent = 0;

% check whether data would is too long for 32-Bit Windows system with 2 GB
% maximum memory, and downsample.
if D.Nchannels*Nblocks*D.Radc > 2^27
    f = ceil(log2(D.Nchannels*Nblocks*D.Radc) - 27); % downsampling factor 2^f
    warning('Data set will be downsampled from %d to %d Hz', D.Radc, D.Radc/(2^f));
    D.Radc = D.Radc/(2^f);
else
    f = 0;
end

% progress bar
spm_progress_bar('Init', Hdata.Head.NRec, 'Seconds converted'); drawnow;
if Hdata.Head.NRec > 100, Ibar = floor(linspace(1, Hdata.Head.NRec,100));
else, Ibar = [1:Hdata.Head.NRec]; end

for t = 1:Nblocks
    
    data = readbdf(Hdata,t,0);
    
    % calibrate the data
    % data.Record = data.Record.*repmat(data.Head.Cal, 1, size(data.Record, 2));

    
    % store data
    %-----------------------------------------------------
    % the EEG data
    d = data.Record(D.channels.eeg, :);
        
    % horizontal EOG
    if Cheog ~= 0
        if length(Cheog) == 2
            % take difference
            d = [d; data.Record(Cheog(1), :) - data.Record(Cheog(2), :)];
        else
            d = [d; data.Record(Cheog(1), :)];
        end        
    end
    
    % vertical EOG
    if Cveog ~= 0
        if length(Cveog) == 2
            % take difference
            d = [d; data.Record(Cveog(1), :) - data.Record(Cveog(2), :)];
        else
            d = [d; data.Record(Cveog(1), :)];
        end              
    end
    
    % reference
    if Creference ~= 0
        if length(Creference) > 1
            % take average
            ref = mean(data.Record(Creference, :));
        else
            ref =  data.Record(Creference(1), :);
        end    
        d = d - repmat(ref, D.Nchannels-no_exg, 1);
    end


    if more_exg ==1
        for kl=1:no_exg
     
          d = [d; data.Record(Chexg(kl), :)];
      end
    end



    % downsampling by simple decimating
    d = d(:, 1:2^f:end);
    fwrite(fpout, d, 'float32');
    Nsamples = size(data.Record, 2); % the original size is needed for events below
    D.Nsamples = D.Nsamples + size(d, 2); % save downsampled Nsamples

    % identify events
    Ezero = 0;
    if ~isempty(Echannel)
        diffs=diff(data.Record(Echannel,:));
        idiffs=find(diffs~=0);
		idiffs=idiffs+1;
        
        % safe-guard against some unexplained zeroing of Echannel
        if data.Record(Echannel, 1) == 0
            warning('Event channel is all zero.');
            Ezero = 1;
            break;
        end
		bytes=dec2bin(data.Record(Echannel,1));
		bytes=bytes(end-7:end);
		bytes=flipdim(bytes,2);
		test_event=bin2dec(bytes);
		if test_event~=lastevent
			idiffs=[1,idiffs];
		end
		idsf=find(diff(idiffs)==1);
		if ~isempty(idsf)
			idiffs(idsf)='';
        end
        

        if ~isempty(idiffs)
			if idiffs(1)==1 & lidiffs==1
				D.events.code=D.events.code(1,1:end-1);
				D.events.time=D.events.time(1,1:end-1);
			end
            for i = idiffs
                bytes=dec2bin(data.Record(Echannel,i));
                bytes=bytes(end-7:end);
                bytes=flipdim(bytes,2);
                event=bin2dec(bytes);
			
                D.events.code = [D.events.code event-128];
                if i==512
					lidiffs=1;
				else
					lidiffs=0;
				end
				lastevent=event;
            end
            D.events.time = [D.events.time idiffs+(t-1)*Nsamples];
        end
    end   
    
    if ismember(t, Ibar)
        spm_progress_bar('Set', t);
        drawnow;
    end

end

if f ~= 0 
   % if downsampling
    D.events.time = ceil(D.events.time/2^f);
end

D.Nevents = 1;
D.events.types = unique(D.events.code);
D.datatype = 'float32';

D.modality = 'EEG';
D.units = '\muV';

fclose(fpout);

D.fname = [F '.mat'];

if spm_matlab_version_chk('7') >= 0
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

spm_progress_bar('Clear');

spm('Pointer','Arrow');
