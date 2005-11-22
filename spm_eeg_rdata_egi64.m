function D = spm_eeg_rdata_egi64(S)
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
% $Id: spm_eeg_rdata_egi64.m 304 2005-11-22 19:43:44Z stefan $

try
	Fdata = S.Fdata;
catch
	Fdata = spm_select(1, '\.txt$', 'Select egitxt-file');
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

Fout = [F '.dat'];
D.fnamedat = ['e_' Fout];
fpd = fopen(fullfile(P, D.fnamedat), 'w');
% read header

fh=fopen([Fdata],'r');
if fh==-1
	error('wrong filename')
end
%read eeg txt file plus labels
for junk=1:3
	line=fgetl(fh);
end
D.events.start=fgetl(fh);;
D.events.start=str2num( D.events.start(7:end))-1;
D.events.stop=fgetl(fh);;
D.events.stop=str2num(D.events.stop(5:end))-1;
D.Nsamples=D.events.stop-D.events.start+1;
fgetl(fh);
D.events.Ntypes=fgetl(fh);
D.events.Ntypes=str2num(D.events.Ntypes(end));
D.Nchannels=fgetl(fh);
D.Nchannels=str2num(D.Nchannels(end-1:end))+1;

%%%% somewhere here must be the sampling rate. Add to D.radc

n=0;
D.events.types=[];
D.events.code=[];
for t=1:D.events.Ntypes
	junk=fgetl(fh);
	D.events.types=[D.events.types,str2num(junk(end))];
	junk=fgetl(fh);
	n=str2num(junk(end-1:end));
	D.events.code=[D.events.code,ones(1,n)*D.events.types(t)];
	junk=fgetl(fh);
	D.Radc=fgetl(fh);
	D.Radc=str2num(D.Radc(end-3:end));
end
try 
    D.events.start = S.events.start;
catch
    D.events.start =...
        spm_input('start of pres-stimulus time [ms]', '+1', 'r', '', 1);
end
D.events.start = ceil(-D.events.start*D.Radc/1000);
D.events.stop = D.Nsamples-D.events.start-1;
D.Nevents=length(D.events.code);
D.channels.Bad=[];
for n=1:D.Nchannels
	D.channels.name{n}=num2str(n);

end
D.channels.eeg=1:D.Nchannels;
D.channels.order=1:D.Nchannels;
D.channels.heog=[];
D.channels.veog=[];

D.channels.reference = 0;
D.channels.ref_name = 'NIL';
D.channels.reference = 65;
D.channels.ref_name = 'Cz';
%%% need to check
D.datatype= 'int16';

fgetl(fh);
fgetl(fh);

gain=fgetl(fh);
gain=str2num(gain(5:end));

fgetl(fh);
fgetl(fh);

zero_no=fgetl(fh);
zero_no=str2num(zero_no(5:end));
fgetl(fh);
fgetl(fh);
fgetl(fh);
spm_progress_bar('Init', (D.Nevents), 'Events read'); drawnow;
if (D.Nevents) > 100, Ibar = floor(linspace(1, (D.Nevents),100));
else, Ibar = [1:D.Nevents]; end
D.scale.dim = [1 3];
D.events.time=[];
for n=1:D.Nevents
	d=zeros(D.Nsamples,D.Nchannels);
	if ismember(n, Ibar)
		spm_progress_bar('Set', n);
		drawnow;
	end
	for time_pt=1:D.Nsamples
		tempdata=fgetl(fh);
		tempdata=str2num(tempdata(13:end));
		tempdata=(tempdata-zero_no).*(400./gain); %corrects for offset and gain
		tempdata(1,65)=0;
		d(time_pt,:)=tempdata;
		
		
	end
	D.events.time=[D.events.time,(n-1)*D.Nsamples+1];
	D.scale.values(:, n) = spm_eeg_write(fpd, d', 2, D.datatype);
	fgetl(fh);
end

D.events.reject=zeros(1,D.Nevents);
spm_progress_bar('Clear');
fclose(fpd);

D.modality = 'EEG';
D.units = '\muV';

D.fname = [F '.mat'];
D.fname = ['e_' D.fname];
if str2num(version('-release'))>=14
	save(fullfile(P, D.fname), '-V6', 'D');
else
	save(fullfile(P, D.fname), 'D');
end

spm_progress_bar('Clear');

spm('Pointer','Arrow');
