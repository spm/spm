function D = spm_eeg_rdata_FIF(S)
% function to read in NeuroMag *.FIF data to SPM5
% FORMAT Do = spm_eeg_rdata_FIF(S)
% 
% S		    - struct (optional)
% (optional) fields of S:
% F			- continuous or averaged FIF file to read
% path      - new path for output
% Pout      - new filename for output
% Fchannels - channel template file
% sclchan   - [0|1] turn channel scaling on or off
% sclfact   - channel scale factors
% HPIfile   - .fif or .mat file (with full path) from which to load HPI points
% conds     - which sets to load from an averaged data file
% twin      - time window (ms) to read, relative to start of continuous
%             data, or event of averaged data
% trig_chan - trigger channel for continuous data
%
%
% NOTES: 
%    - Requires the Fiff Access toolbox written by Kimmo Uutela
%      http://www.kolumbus.fi/~w132276/programs/meg-pd/
%    - Fiff Access does not handle long int, so may need to convert 
%      to short int or float first (eg using MaxFilter)
%    - Does NOT apply any SSP vectors within FIF file
%    - Requires separate function spm_eeg_rdata_FIF_channels.m
%      This function returns channel information in Brainstorm format,
%      only part of which is stored in SPM's D structure.
%
% Will put SPM5 *.mat file in same directory as FIF file unless 
% alternative S.path or S.Pout passed
%
% Rik Henson (5/06/07), with thanks to Danny Mitchell and Jason Taylor

try
    	Fdata = S.Fdata;
catch
    	Fdata = spm_select(1, '\.fif$', 'Select FIF file');
end
P = spm_str_manip(Fdata, 'H');

try
    	D.path = S.path;
catch
    	D.path = P;
end

if exist(Fdata)==2
    if loadfif(Fdata,'sets')==0
	rawflag = 1;
    else
	rawflag = 0;
    end
else
    error('Failed to find %s',Fdata);
end

try
    	Fchannels = S.Fchannels;
catch
    	Fchannels = spm_select(1, '\.mat$', 'Select channel template file', {}, fullfile(spm('dir'), 'EEGtemplates'));
end
D.channels.ctf = spm_str_manip(Fchannels, 't');

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','MEG data conversion ',0);

% Could read bad channels using "badchans", if could get "megmodel" to work!
%[bad,names] = badchans;

D.channels.Bad = [];

% compatibility with some preprocessing functions
D.channels.heog = 0;
D.channels.veog = 0;
D.channels.reference = 0;
D.channels.ref_name = 'none';

%%%%%%%%%%%%%% Read channel information

[Channel,B] = spm_eeg_rdata_FIF_channels(Fdata);
D.Nchannels = length(Channel);

% Reorder channels by channel type (magnotometer vs gradiometer) 
% so that subsequent display in SPM is easier to select/deselect!

senstype = strvcat(Channel(:).SensType);
mags = strmatch('Magnetometer',senstype);
grds = strmatch('Gradiometer',senstype);
reord = [mags; grds];
Channel = Channel(reord);

D.channels.eeg = 1:D.Nchannels;
D.channels.name = cat(1,Channel(:).Name);

Weights = cat(2,Channel(:).Weight)';
D.channels.Weight = reshape(Weights,2,length(Weights)/2)';

% Scale the data by the weights, unless the user overrides 
% with S.sclchan = 0 or with only scale factor in S.sclfact

try 
     sclchan = S.sclchan;
catch
     sclchan = 1;
end

if sclchan
    try
	sf = S.sclfact(:);
        if length(sf) ~= D.Nchannels
	    error(sprintf('Only %d scale factors supplied; %d expected',length(sf),D.Nchannels))
	else
	    disp('Reordering scale factors by mags then grds')
	    D.channels.scaled = sf(reord);
	end
     catch
	disp('Using weights from channel file to scale data')
     	D.channels.scaled = D.channels.Weight(:,1);
	D.channels.Weight = D.channels.Weight./repmat(D.channels.scaled,1,size(D.channels.Weight,2));
     end
     disp('Will rescale channel data by...')
     disp(D.channels.scaled')
else
     D.channels.scaled = ones(D.Nchannels,1);
end

% Sensor positions and orientations not normally in D structure,
% but perhaps they should be!!! 
% (NB: These are properties of sensor array, ie helmet, so do not 
% depend on the subject, but their coordinates are in DEVICE space, 
% and need to be rotated to HEAD space, which does depend on subject)

Tmat = loadtrans(Fdata,'DEVICE','HEAD');

% Extract from Brainstorm 6xN format (using Brainstorm fields)....
CoilLoc = cat(2,Channel(:).Loc)';
CoilOrt = cat(2,Channel(:).Orient)';
GradOrt = cat(2,Channel(:).GradOrient)';

% Transform first coil info into head coordinates, and convert to mm
Loc = CoilLoc(:,1:3);
Ort = CoilOrt(:,1:3);
Loc = Tmat * [Loc ones(size(Loc,1),1)]';
Ort = Tmat * [Ort zeros(size(Ort,1),1)]';
D.channels.Loc = 1000*Loc(1:3,:);
D.channels.Orient = Ort(1:3,:);

% Use first coil for location and orientation of sensor itself (this is correct for mags; slight displacement for planar grads, but this only for display purposes; full coil location and orientation is retained in Loc and Orient fields)
D.channels.pos3D = D.channels.Loc(1:3,:)';
D.channels.ort3D = Tmat * [GradOrt(:,1:3) zeros(size(GradOrt,1),1)]';
D.channels.ort3D = D.channels.ort3D(1:3,:)';

% Transform second coil info into head coordinates (second coil in mags is NaN)
Loc = CoilLoc(:,4:6);
Ort = CoilOrt(:,4:6);
Loc = Tmat * [Loc ones(size(Loc,1),1)]';
Ort = Tmat * [Ort zeros(size(Ort,1),1)]';
D.channels.Loc = [D.channels.Loc; 1000*Loc(1:3,:)];
D.channels.Orient = [D.channels.Orient; Ort(1:3,:)];

% Could add option to get sensor locations from template, if not found in file
%D.channels.pos3D = Tmat * [Csetup.SensLoc ones(size(Csetup.SensLoc,1),1)]';
%D.channels.ort3D = Tmat * [Csetup.SensOr zeros(size(Csetup.SensOr,1),1)]';


% Read in Head Position Points (HPI) from Isotrak in HEAD space
% (If MaxMove has been used, these may have been deleted, so allow them to 
% be loaded seperately from a reference fif or mat file provided in S.HPIfile. 
% djm 27/6/07)
co=[];
try 
    [pth name ext]=fileparts(S.HPIfile);
    if ext=='.fif'
        [co,ki,nu] = hpipoints(S.HPIfile);
    else
        load(S.HPIfile,'co','ki','nu');
    end
catch
    [co,ki,nu] = hpipoints(Fdata);
end
if isempty(co) | ~exist('ki') | ~exist('nu') | length(co)<7
    error('Failed to load HPI points!')
end

D.channels.fid_eeg = 1000*co(:,find(ki==1))';
[dummy,nasion] = max(D.channels.fid_eeg(:,2));
[dummy,LE] = min(D.channels.fid_eeg(:,1));
[dummy,RE] = max(D.channels.fid_eeg(:,1));
D.channels.fid_eeg = D.channels.fid_eeg([nasion LE RE],:);
disp('Reordering cardinal points to...'); disp(D.channels.fid_eeg);
D.channels.fid_coils = 1000*co(:,find(ki==2))';
D.channels.headshape = 1000*co(:,find(ki==4))';

try 
	pflag = S.pflag
catch
	pflag = 1;
end
if pflag
 xyz = D.channels.pos3D;
 ori = D.channels.ort3D;
 fid = D.channels.fid_eeg;
 hsp = D.channels.headshape;
 nam = strvcat(D.channels.name{1:D.Nchannels});
 h    = spm_figure('GetWin','Graphics');
 clf(h); figure(h), hold on
 lstr = xyz(:,1:3)';
 lend = lstr+ori(:,1:3)'*10;
 plot3(xyz(:,1),xyz(:,2),xyz(:,3),'o');
% t=text(xyz(1:102,1),xyz(1:102,2),xyz(1:102,3),nam(1:102,4:7));
 t=text(xyz(:,1),xyz(:,2),xyz(:,3),nam(:,4:7));
 set(t,'FontSize',6);
 line([lstr(1,:); lend(1,:)],[lstr(2,:); lend(2,:)],[lstr(3,:); lend(3,:)]);
 plot3(fid(:,1),fid(:,2),fid(:,3),'g+')
 plot3(hsp(:,1),hsp(:,2),hsp(:,3),'r.')
 title('Planar Grad locations slightly displaced for visualisation')
end
rotate3d on

%%%%%%%%%%%%%% Find channels in template for display

Csetup = load(D.channels.ctf);

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


if rawflag == 0
%-----------------------------------
% IMPORT AVERAGE DATA:

 [B.totconds,B.comments] = loadfif(Fdata,'sets');
 disp(B.comments)

 try
	conds = S.conds;
 catch
	conds = spm_input(sprintf('Which conditions? (1-%d)',B.totconds),'+1','r',[1:B.totconds]);
 end

 D.Nevents = length(conds);
 D.events.types = [1:D.Nevents];
 D.events.Ntypes = length(D.events.types);
 D.events.code = [1:D.Nevents];
 D.events.reject = zeros(1, D.events.Ntypes);
 D.events.repl = ones(1, D.events.Ntypes);	% (lost original number of trials?)

%%%%%%%%%%%%%% Read data

 for c = conds			% Currently assumes B.t0 same for all conds
  [B.data{c},B.sfreq,B.t0] = loadfif(Fdata,c-1,'any');
  elen(c) = size(B.data{c},2);
 end
 D.Radc = B.sfreq;
 elen(elen==0)=[]; % djm 27/06/07
 disp('epoch lengths (in samples)...'),disp(elen)

%%%%%%%%%%%%%% Define epoch (S.twin, if specified, is eg [-100 300])

 B.t0 = round(B.t0*1000);	% convert to ms  % added round. djm 27/6/07. 
 try	
	twin = S.twin;
 catch
	swin = [1 min(elen)];
	twin = round((swin-1)*1000/D.Radc + B.t0);
	if any(diff(elen))
	  twin = spm_input('Sample window? (ms)','+1','r',twin,2,twin)';
	end
 end
	
 swin = round((twin-B.t0)*D.Radc/1000)+1;
 if length(twin)~=2 | twin(1) < B.t0 |  twin(1) > -1 | swin(2) > min(elen)
	error('twin outside range for all conditions')
 end

 D.Nsamples = swin(2) - swin(1) + 1;
 D.events.start = -round(twin(1)*D.Radc/1000);
 D.events.stop = swin(2) - swin(1) - D.events.start;

%%%%%%%%%%%%%% Reformat data
% !!Could do baseline correction here, particularly if baseline period changed

 for c = 1:length(conds)
   d(:,:,c) = B.data{conds(c)}(B.chanfilt,swin(1):swin(2));
 end

 d = d(reord,:,:);

 d = d.*repmat(D.channels.scaled,[1 size(d,2) size(d,3)]);

try
    % option to provide different output file in S.Pout - djm 27/6/07
    [f1, f2, f3] = fileparts(S.Pout);
    D.fname = [f2 '.mat'];      
    D.fnamedat = [f2 '.dat']; 
    if ~isempty(f1); D.path = f1; end;
catch
    [dummy,stem,ext] = fileparts(Fdata);
    D.fname = strcat('me_',stem,'.mat');
    D.fnamedat = strcat('me_',stem,'.dat');
end

else
%-----------------------------------
% IMPORT RAW DATA:

%%%%%%%%%%%%%% Find trigger channel

 D.Nevents = 1;
 try 
	trig_chan = S.trig_chan;
 catch
	trig_chan = spm_input('Trigger channel name?(0=none)','+1','s','STI101');
 end
 disp(sprintf('Trigger channel %s',trig_chan))

 if ~strcmp(trig_chan,'0')
  k = [];
  while isempty(k)
   c = 1;
   while isempty(k) & c <= size(B.chans,1)
	k = findstr(B.chans{c},trig_chan);
	c = c+1;
   end
   if c > size(B.chans,1)
	disp(B.chans)
	disp('Channel not found...')
	trig_chan = spm_input('Trigger channel name?(0=none)?','+1','s','STI101');
   end
  end

%%%%%%%%%%%%%% Read events

  spm('Pointer', 'Watch'); drawnow;
  [B.data, D.Radc] = rawchannels(Fdata,trig_chan);

  D.Nsamples = size(B.data,2);
% If assume trigger onset always positive...
%  pos = find(B.data>0); trig = find(diff(pos)>1); D.events.time = pos(trig);
  trig = find(abs(diff(B.data))>0)+1;
  D.events.time = trig(1:2:end);	% exclude offsets (assumes duration >1 sample)
  disp(sprintf('%d triggers found...',length(D.events.time)))
  D.events.code = B.data(D.events.time);
  D.events.types = unique(D.events.code);
  disp(sprintf('Trigger codes = %s',mat2str(D.events.types)))
  D.events.Ntypes = length(D.events.types);

 else 
% inefficient to read one channel to get Nsamples, but hey....!
  spm('Pointer', 'Watch'); drawnow;
  [B.data, D.Radc] = rawchannels(Fdata,B.chans{1});
  D.Nsamples = size(B.data,2);
  D.events.time = 1;
  D.events.code = 1;
  D.Nsamples = 0;	% because don't know length yet
 end

 spm('Pointer', 'Arrow'); drawnow;

%%%%%%%%%%%%%% Select sample
 try	
	twin = S.twin;
 catch
	swin = [1 D.Nsamples];
	twin = (swin-1)/D.Radc;
	twin = spm_input('Sample window? (s)','+1','r',twin,2,twin)';
 end

% If want to prespecify all, without precise knowledge of length
 if twin(2) == Inf	
	swin(1) = round(twin(1)*D.Radc)+1;
	swin(2) = D.Nsamples;
 else
	swin = round(twin*D.Radc)+1;
 end
 
 if length(twin)~=2 | twin(1) < 0 | swin(2) > D.Nsamples
	error('sample outside range')
 end

%%%%%%%%%%%%%% Read raw data

 spm('Pointer', 'Watch'); drawnow;

 spm_progress_bar('Init', 100, 'Samples read'); drawnow;
 Ibar = floor(linspace(1, swin(2), 100));
 status=1; lfi=0;
 
 B.data=[];
 rawdata('any',Fdata);
 while status
	[b,status] = rawdata('next');
	B.data = [B.data b];
	fi = find(Ibar-size(B.data,2) > 0);
        if isempty(fi)
	  status = 0;
	elseif fi(1)>lfi
	  lfi=fi(1);
          spm_progress_bar('Set', lfi);
          drawnow;
        end
 end
 rawdata('close');

 spm_progress_bar('Clear');

 D.Nsamples = swin(2) - swin(1) + 1;
 disp(sprintf('%d samples read (%f seconds)',D.Nsamples,D.Nsamples/D.Radc))

%%%%%%%%%%%%%% Reformat data

 d = B.data(B.chanfilt,swin(1):swin(2));

 d = d(reord,:);

 d = d.*repmat(D.channels.scaled,1,size(d,2));

 try
    % option to provide different output file in S.Pout - djm 27/6/07
    [f1, f2, f3] = fileparts(S.Pout);
    D.fname = [f2 '.mat'];      
    D.fnamedat = [f2 '.dat'];      
    if ~isempty(f1); D.path = f1; end;
 catch
    [dummy,stem,ext] = fileparts(Fdata);
    D.fname = strcat(stem,'.mat');
    D.fnamedat = strcat(stem,'.dat');
 end
 
end

d = d*10^15;
D.units = 'fT';
D.modality = 'MEG';


%%%%%%%%%%%%%% Any EOG?
% Hack: assumes any EEG channels are the two EOG channels!!!
% Hack: if there are two EEG (EOG) channels, it assumes VEOG first!!

B.eogfilt = find(B.chtypes==2 | B.chtypes==202)

if ~isempty(B.eogfilt)

  if length(B.eogfilt)~=2
   error('Discovered other than 0 or 2 EEG channels???')
  end

  disp(sprintf('Found two EEG channels %s and %s - assuming are VEOG and HEOG respectively!!!',B.chans{B.eogfilt(1)},B.chans{B.eogfilt(2)}))

  if rawflag == 0
   for c = 1:length(conds)
    d(D.Nchannels+[1:2],:,c) = B.data{conds(c)}(B.eogfilt,swin(1):swin(2))*10^6;
   end
  % postmultiplier is to convert to uV
  else
    d(D.Nchannels+[1:2],:) = B.data(B.eogfilt,swin(1):swin(2))*10^6;
  end

  D.channels.veog = D.Nchannels+1;
  D.channels.name{D.channels.veog} = 'VEOG';
  D.channels.heog = D.Nchannels+2;
  D.channels.name{D.channels.heog} = 'HEOG';
  D.Nchannels = D.Nchannels+2;

  for i = [D.channels.veog D.channels.heog]
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
end



%%%%%%%%%%%%%% Write data

D.scale = ones(D.Nchannels, 1, D.Nevents);
D.datatype  = 'float32';

fpd = fopen(fullfile(D.path, D.fnamedat), 'w');

for e = 1:D.Nevents
%    dd = squeeze(d(:, :, i));
%    D.scale(:, 1, i) = spm_eeg_write(fpd, dd, 2, D.datatype);
    for s = 1:D.Nsamples	
	fwrite(fpd, d(:,s,e), 'float');
    end
end

fclose(fpd);

if str2num(version('-release'))>=14 
    save(fullfile(D.path, D.fname), '-V6', 'D');
else
    save(fullfile(D.path, D.fname), 'D');
end

spm('Pointer', 'Arrow');




