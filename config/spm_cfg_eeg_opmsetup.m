function setup = spm_cfg_eeg_opmsetup
% Configuration file for M/EEG OPM set up
%__________________________________________________________________________
% Copyright (C) 2017-2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney, Gareth Barnes
% $Id: spm_cfg_eeg_opmsetup.m 7429 2018-09-28 09:29:20Z spm $


R        = cfg_files;
R.tag    = 'R';
R.name   = 'OPM to ADC file (specific to OPM connections to hardware)';
R.filter = 'txt';
R.num    = [1 1];
R.help   = {'Select the recording set up file- this file contains a list of OPM channel names and which binary data (eg labview) channels they were linked to'};

C        = cfg_files;
C.tag    = 'C';
C.name   = 'Scanner-cast file (specific to specific subject and scanner-cast)';
C.filter = 'txt';
C.num    = [1 1];
C.help   = {'Select the scanner cast defining file- this file contains positions and orientations of OPMs in native MRI coords, along with a list of names for the slots'};

RC        = cfg_files;
RC.tag    = 'RC';
RC.name   = 'OPM to cast file';
RC.filter = 'txt';
RC.num    = [1 1];
RC.help   = {'Select the file which shows which OPM channels were in which scanner-cast slots'};



M        = cfg_files;
M.tag    = 'M';
M.name   = 'Native space MRI of subject';
M.filter = {'nii','img'};
M.num    = [1 1];
M.help   = {'The native space MRI for this subject'};

fid         = cfg_entry;
fid.tag     = 'fid';
fid.name    = 'Fiducial points in MRI [nas;lpa;rpa] in mm';
fid.strtype = 'r';
fid.num     = [3 3];
fid.help    = {'Enter approx fiducial points from native space mri in format [nas;lpa;rpa] (mm) (these will be saved in a BIDs format file for later)'};



L        = cfg_files;
L.tag    = 'lbv';
L.name   = 'Binary file containing raw OPM data';
L.filter = {'lvm','dat'};
L.num    = [1 1];
L.help   = {'Select the labview or binary data file containing the raw data. '};

binind         = cfg_entry;
binind.tag     = 'binind';
binind.name    = 'Binary data channel indices';
binind.strtype = 'r';
binind.num     = [1 Inf];
binind.help    = {'Data channel indices in binary file'};



% Trigger channels
convert         = cfg_entry;
convert.tag     = 'convert';
convert.name    = 'Conversion from ADC units to Tesla';
convert.strtype = 'r';
convert.num     = [1 1];
convert.help    = {'Conversion from ADC units to T (typically 1e6/2.7)'};
convert.val=    {1e6/2.7};

% Trigger channels
trigs         = cfg_entry;
trigs.tag     = 'trigind';
trigs.name    = 'Trigger channel indices';
trigs.strtype = 'r';
trigs.num     = [1 Inf];
trigs.help    = {'Trigger channel indices in binary file'};


prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Output prefix';
prefix.strtype = 's';
prefix.help    = {'prefix for SPM object'};
prefix.val= {'OPM'};

band         = cfg_entry;
band.tag     = 'band';
band.name    = 'Filter unepoched data to';
band.strtype = 'r';
band.num     = [2 1];
band.help    = {'Frequency band in Hz (2nd order butterworth, 2 way)'};

epwin         = cfg_entry;
epwin.tag     = 'epwin';
epwin.name    = 'Time window around triggers (ms)';
epwin.strtype = 'r';
epwin.num     = [2 1];
epwin.help    = {'Epoch window relative to triggers'};

ncan         = cfg_entry;
ncan.tag     = 'ncan';
ncan.name    = 'reference noise cancellation';
ncan.strtype = 'r';
ncan.num     = [2 1];
ncan.val  = { [ 1 1] };
ncan.help    = {'Flags (0 or 1) for reference noise cancellation. [1 0] use reference derivatives; [ 0 1] use global signal; [1 1] use both'};


setup          = cfg_exbranch;
setup.tag      = 'OPMsetup';
setup.name     = 'Set up OPM recording';
setup.val      = {R, C, RC, M, fid, L, convert, trigs,band,epwin, ncan, prefix};
setup.help     = {'Set up OPM data in spm format'};
setup.prog     = @opmsetup;
setup.vout     = @vout_opmsetup;


%==========================================================================
%-function out = opmsetup(job)
%==========================================================================
function out = opmsetup(job)

%-READ IN RECORDING SET UP FILE 
%==========================================================================
% this table has defines which OPM sensor names go to which ADC channels
% format: adchan\t opmname\n
% use some default like XX for no channel connected, or delete redundant
% rows
recordingsetup=cell2mat(job.R);
fprintf('\n Reading basic recording set up file %s', recordingsetup);
matLabel = readtable(recordingsetup,'ReadVariableNames',0,'Delimiter','\t');  %% This links the recorded data to sensor names
adcind=matLabel{:,1};
fprintf('\nFound a total of %d listed connections\n',length(adcind));
adcOPMs=table2array(matLabel(:,2)); % OPM names attached to the ADC channels
timeind=strmatch('Time',adcOPMs,'exact');
if isempty(timeind)
    warning('No time channel found, assuming it is channel 1');
else
    fprintf('\n Found time channel on %d',timeind);
end

%-READ IN SCANNER-CAST SPECIFYING FILE
%==========================================================================
% gives where the slots are in native space and which OPM channels were in these slots
% each slot has position and orienttation columns 1:3 columns 4:6 plus a column 7
castsetup=cell2mat(job.C);
fprintf('\n Reading scanner cast file %s ', castsetup);
cast = readtable(castsetup,'ReadVariableNames',0,'Delimiter','\t'); %% this is a scanner-cast specific file
slotname=table2array(cast(:,7));
slotname=strtrim(slotname);
Ncastslots=size(cast,1); %% number of possible channels in cast
posor=table2array(cast(:,1:6));
fprintf('\nFound %d OPM slots on headcast\n',Ncastslots);
fprintf('\n Names ');


%-NOW READ IN OPM2SLOT FILE
%==========================================================================
% this is the file that will change most often
% and links the opms to specific slots in scanner-cast
% format : opmname\t slotname\n
% use REF in 2nc column to define a reference channel (i.e. used but not in cast)
opm2slotfile=cell2mat(job.RC);
opm2slot = readtable(opm2slotfile,'ReadVariableNames',0,'Delimiter','\t'); %% this how the OPM channels were slotted into the cast
opm2slot=table2array(opm2slot);
fprintf('\n*****************************\n')
fprintf('\nOPM\t\tSlot\tADC chan\n')
for f=1:size(opm2slot,1)
    fprintf('%s\t\t',opm2slot{f,1});
    
    if ~(strcmp(opm2slot{f,2},'REF'))
        slotind(f)=strmatch(opm2slot{f,2},slotname); %% find the OPM slot name in the scannercastfile (with positions and orientations)
        if isempty(slotind)
            error('Could not find slot %s in cast file',opm2slot(f,2));
        end
        fprintf('%s\t\t',slotname{slotind(f)});
    else
        slotind(f)=-1; %% no scanner-cast slot as this is a reference
        fprintf('REF\t\t');
    end
    listind=strmatch(opm2slot{f,1},adcOPMs,'exact'); %% find OPM channel in list of ADC connections (in pinout file)
    if isempty(listind)
        error('Could not find OPM %s in OPM to ADC file',opm2slot(f,1));
    end
    opmadcind(f)=adcind(listind); %% adc channel corresponding to opm opm2slot(f,1)
    fprintf('%d\n',opmadcind(f));
end
fprintf('\n*****************************\n')
megind=find(slotind>0); %% indices of OPM chans in headcast


%-READ IN LABVIEW BINARY FILE
%==========================================================================
% containing the actual recorded data for OPMs and trigger channels
labviewfile=cell2mat(job.lbv);
[a1,b1,c1]=fileparts(labviewfile);
fprintf('\nReading binary file %s', labviewfile);
[time,B]=spm_opm_read_lvm(labviewfile,timeind);
srate=1./mean(diff(time));
fprintf('\nRead %d samples of %d channels at %3.2fHz sampling rate\n',size(B,1),size(B,2),srate);
% every second column is a radial sensor

alltrigs=round(job.trigind); %% the adc channels containing the triggers

megData = B(:,[opmadcind alltrigs])'; %% get data from the adcchannels that the OPMs were connected to.


outbase=sprintf('spm_%s_%s.dat',job.prefix,b1); %% name of output file


%-Convert to SPM
%--------------------------------------------------------------------------
if exist(outbase)
    warning('SPM file already exists, OVERWRITING');
    D=spm_eeg_load(outbase);
    D.delete;
end

D = spm_opm_convert(megData,outbase,srate,job.convert);


%-Set appropriate Channel labels and types
%--------------------------------------------------------------------------
Nopmchans=size(opm2slot,1);
D = chanlabels(D,1:Nopmchans,opm2slot(:,1));
D = chantype(D,megind,'MEG');
D = chantype(D,find(slotind==-1),'Ref');
D=chantype(D,length(opmadcind)+1:length(opmadcind)+length(alltrigs),'Trig')
for f=1:length(alltrigs)
    D=chanlabels(D,length(opmadcind)+1:length(opmadcind)+length(alltrigs),sprintf('tr%d',alltrigs(f)));
end
save(D);

%-Slot positions and orientation into MEG object
%==========================================================================
% this table maps positions/orientations to headcast numbers and sensor
% labels

if length(megind)~=size(posor,1)
    error('File size mismatch');
end
figure; hold on;
grad.label = opm2slot(megind,1); %% names of the MEG channels in headcast
megpos=[];meglabels=[];
for f=1:length(megind) %% megind indexes into slotind which gives the index into posor
    castind=slotind(megind(f));
    grad.coilpos(f,:) = posor(castind,1:3);
    grad.coilori(f,:) = posor(castind,4:6);
    %plot3(grad.coilpos(f,1),grad.coilpos(f,2),grad.coilpos(f,3),'r.');
    %text(grad.coilpos(f,1),grad.coilpos(f,2),grad.coilpos(f,3),grad.label{f});
end


grad.tra = eye(numel(grad.label));



grad.chanunit = repmat({'T'}, numel(grad.label), 1); %% grad units always in Tesla (data units in fT)
grad = ft_datatype_sens(grad, 'amplitude', 'T', 'distance', 'mm');

D = sensors(D, 'MEG', grad);

%-ORIENT SENSORS IN 2D BASED ON MEAN ORIENTATION
%==========================================================================
n1=mean(grad.coilori); n1= n1./sqrt(dot(n1,n1));
t1=cross(n1,[0 0 1]);
t2=cross(t1,n1);

for f=1:length(megind)
    pos2d(f,1)=dot(grad.coilpos(f,:),t1);
    pos2d(f,2)=dot(grad.coilpos(f,:),t2);
    plot(pos2d(f,1),pos2d(f,2),'r.');
    text(pos2d(f,1),pos2d(f,2),grad.label{f});
end

S=[]; %% SET SENSOR 2D POS
S.D=D;
S.xy= pos2d';
S.label=grad.label;
S.task='setcoor2d';
D=spm_eeg_prep(S);
D.save;
pos=D.coor2D';
% for f=1:length(megind),
%     
%     h=text(pos(f,1),pos(f,2)+10,grad.label{f});
%     set(h,'color','m');
% end;

D.save;

%-NOW COREGISTER
%==========================================================================
% easiest to do this here as the opms and head are defined in same native mri space
% working from native to native space so fiducials should give identity transform

newfid=[];
newfid.fid.label = {'nas', 'lpa', 'rpa'}';
newfid.fid.pnt=job.fid;

mri_fids=newfid.fid.pnt; %%% the MRI fiducials will be the same as the MEG as everything is defined in native space for OPMs
D = fiducials(D, newfid);
save(D);

%-now write out a bids format file to be read by coreg procedure
%--------------------------------------------------------------------------

[a1,b1,c1]=fileparts(cell2mat(job.M));
json_fid = [];

%json_fid.IntendedFor = cell2mat(fullfile(job.M));
json_fid.IntendedFor = [b1 c1]; %% does not like full path for some reason
json_fid.AnatomicalMRICoordinateSystem = 'other';
json_fid.AnatomicalMRICoordinateUnits  = 'mm';
json_fid.CoilCoordinates.nas = spm_vec(mri_fids(1, :));
json_fid.CoilCoordinates.lpa    = spm_vec(mri_fids(2, :));
json_fid.CoilCoordinates.rpa    = spm_vec(mri_fids(3, :));
json_fid.CoilCoordinateSystem   = 'native';
json_fid.CoilCoordinateUnits    = 'mm';

jsonname=[a1 filesep b1 '_fid.json'];
spm_jsonwrite(jsonname, json_fid);

%test=spm_jsonread(jsonname);



figure;

plot3(posor(:,1),posor(:,2),posor(:,3),'ro'); hold on;
text(posor(:,1),posor(:,2),posor(:,3),slotname);
plot3(mri_fids(:,1),mri_fids(:,2),mri_fids(:,3),'m*');
text(mri_fids(1,1),mri_fids(1,2),mri_fids(1,3),'nasion');
text(mri_fids(2,1),mri_fids(2,2),mri_fids(2,3),'lpa');
text(mri_fids(3,1),mri_fids(3,2),mri_fids(3,3),'rpa');
fprintf('\n Will need to coregister later using the BIDS json file created here')

figure;



%-bandpass filter
%--------------------------------------------------------------------------
fprintf('\n Filtering %d to %d Hz \n', job.band(1), job.band(2));
S = [];
S.D = D;
S.type = 'butterworth';
S.band = 'bandpass';
S.freq = [job.band(1) job.band(2)];
S.dir = 'twopass';
S.order = 2;
nD = spm_eeg_filter(S);


%%%%%%%%%%%%%%%%%epoching %%%%%%%%%%%%%%

% this will be different for different datasets

ev=[];cond={};

figure; % ASSUME ONE BIT PER TRIG CHANNEL FOR NOW AND THRESHOLD AT 2* std deviation

for f=1:length(alltrigs)
    subplot(length(alltrigs),1,f); 
    thresh=std((B(:,alltrigs(f))))*2;
    trigs=find(diff(B(:,alltrigs(f))>thresh)==1)+1;
    plot(D.time,B(:,alltrigs(f))',D.time,ones(size(D.time))*thresh,':',D.time(trigs),ones(size(trigs))*thresh*2,'r*');
    legend(num2str(alltrigs(f)),'thresh','events'); xlabel('time');ylabel('adc value')
    title('trigger');
    for j=1:length(trigs)
        cond{j+length(ev)}=num2str(alltrigs(f));
    end
    ev = [ev; trigs]; %% trigger samples
end

fprintf('\n Epoching');
offsettime=job.epwin(1)/1000;
duration=diff(job.epwin)/1000;
offsetsamples=round(offsettime.*nD.fsample);
durationsamples=round(duration.*nD.fsample);
begsample=ev+offsetsamples;
offset=offsetsamples.*ones(size(begsample));
endsample=begsample+durationsamples;
% should be in format : begin sample, end sample, offset
trl = round([begsample  endsample  offset]);
S = [];
S.D = nD;
S.trl = trl;
S.conditionlabels =cond;
S.bc = 0;
S.prefix = 'e_';
eD = spm_eeg_epochs(S);


%-Reference noise cancellation
%--------------------------------------------------------------------------
refInd = selectchannels(eD,'Ref');
%%spm_opm_denoise(D,refD,derivative,gs,update,prefix)
fprintf('\n Reference noise cancellation with flags %d %d', job.ncan(1),job.ncan(2));
dD =  spm_opm_denoise(eD,eD(refInd,:,:),job.ncan(1),job.ncan(2),1,'d');


save(dD);


out.D = {fullfile(dD.path, dD.fname)};
out.json={jsonname};

fprintf('\n Finished setup');


%==========================================================================
%-function dep = vout_opmsetup(job)
%==========================================================================
function dep = vout_opmsetup(job)
% return dependencies
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Prepared M/EEG Data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});
%dep.tgt_spec   = cfg_findspec({{'strtype','e'}});


dep(2)            = cfg_dep;
dep(2).sname      = 'M/EEG fiducial information';
dep(2).src_output = substruct('.','json');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});

% if job.phase
%     dep(3)            = cfg_dep;
%     dep(3).sname      = 'M/EEG time-frequency phase dataset';
%     dep(3).src_output = substruct('.','Dtph');
%     % this can be entered into any evaluated input
%     dep(3).tgt_spec   = cfg_findspec({{'strtype','e'}});
%
%     dep(4)            = cfg_dep;
%     dep(4).sname      = 'M/EEG time-frequency phase dataset';
%     dep(4).src_output = substruct('.','Dtphname');
%     % this can be entered into any file selector
%     dep(4).tgt_spec   = cfg_findspec({{'filter','mat'}});
% end
%
