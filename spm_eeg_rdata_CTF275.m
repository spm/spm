function D = spm_eeg_rdata_CTF275(S)
%%%% function to read in CTF data to Matlab

try
	timewindow = S.tw;
catch
	timewindow = spm_input('do you want to read in all the data','+1','yes|no',[1 0]);
end
if timewindow ==1 
    timeperiod='all';
else
    try
        timeperiod=S.timeperiod;
    catch
        [Finter,Fgraph,CmdLine] = spm('FnUIsetup','MEG data conversion ',0);
        str = 'time window';
        YPos = -1;
        while 1
            if YPos == -1
                YPos = '+1';
            end
            [timeperiod, YPos] = spm_input(str, YPos, 'r');
            if timeperiod(1) < timeperiod(2), break, end
            str = sprintf('window must increase with time');
        end
    end
end
try
   pre_data = ctf_read(S.Fdata,[],timeperiod,[],0);
catch
 error('wrong folder name')
end

try
    Fchannels = S.Fchannels;
catch
    Fchannels = spm_select(1, '\.mat$', 'Select channel template file', {}, fullfile(spm('dir'), 'EEGtemplates'));
end


D.channels.ctf = spm_str_manip(Fchannels, 't');

D.channels.Bad = [];

% compatibility with some preprocessing functions
D.channels.heog = 0;
D.channels.veog = 0;
D.channels.reference = 0;

I = STRMATCH('UPPT0',pre_data.sensor.label)
if isempty(I)
    warning(sprintf('No parallel port event channel was found in the CTF file: your data will be read without events'))
    D.events.time=[];
    D.events.code=[];
else
    if length(I)==1
        PP1=squeeze(pre_data.data(:,I(1)));
    else
        PP1=squeeze(pre_data.data(:,I(1)));
        PP2=squeeze(pre_data.data(:,I(2)));
    end
    D.events.time=[];
    D.events.code=[];
    inds=find(diff(PP1)>0);
    if ~isempty(inds)
        if length(PP1)<inds(end)+2
            inds(end)='';
        end
        D.events.code=PP1(inds+2)'; %changed to +2 from +1 to avoid errors when changing event code without passing by zero.
        D.events.time=inds'+1;
    end
    if length(I)>1
        inds=find(diff(PP2)<0);
        if ~isempty(inds)
            D.events.code=[D.events.code,PP2(inds+1)'+255];
            
            D.events.time=[D.events.time,inds'+1];
        end
    end
    
    [X,I]=sort(D.events.time);
    D.events.time=D.events.time(I);
    D.events.code=D.events.code(I);
end
sens=strmatch('M',pre_data.sensor.label);
D.channels.name=pre_data.sensor.label(sens);
D.channels.order=[1:length(sens)];
D.Nchannels=length(sens);
D.channels.eeg=[1:length(sens)];
D.Radc=pre_data.setup.sample_rate;
D.Nsamples=pre_data.setup.number_samples;
D.Nevents=pre_data.setup.number_trials;
[pathstr,name,ext,versn]=fileparts(pre_data.folder);
D.datatype= 'float';
D.fname=[name,'.mat'];
D.path=pwd;
D.fnamedat=[name,'.dat'];
if D.Nevents>1
    D.events.start=pre_data.setup.pretrigger_samples;
    D.events.stop=D.Nsamples-pre_data.setup.pretrigger_samples-1;
    D.events.reject=zeros(1,D.Nevents);
    D.events.code=ones(1,D.Nevents);
    D.events.types=1;
    D.events.Ntypes=1;
end
D.scale = ones(D.Nchannels, 1, D.Nevents);

fpd = fopen(fullfile(D.path, D.fnamedat), 'w');
for ev=1:D.Nevents
    for n=1:D.Nsamples
	
		fwrite(fpd, pre_data.data(n,sens,ev).*1e15, 'float');
	
    end
end


    
fclose(fpd);

% --- Save coil/sensor positions and orientations for source reconstruction (in mm) ---

% - channel locations and orientations
SensLoc = [];
SensOr  = [];
for i = 1:length(pre_data.sensor.location);
    if any(pre_data.sensor.location(:,i)) & pre_data.sensor.label{i}(1) == 'M'
        SensLoc = [SensLoc; pre_data.sensor.location(:,i)'];
        SensOr  = [SensOr ; pre_data.sensor.orientation(:,i)'];
    end
end
SensLoc = 10*SensLoc; % convertion from cm to mm
if length(SensLoc) > 275
    warning(sprintf('Found more than 275 channels!\n'));
end

[pth,nam,ext]  = fileparts(D.fname);
fic_sensloc    = fullfile(D.path,[nam '_sensloc.mat']);
fic_sensorient = fullfile(D.path,[nam '_sensorient.mat']);
save(fic_sensloc, 'SensLoc');
save(fic_sensorient, 'SensOr');
clear SensLoc

% for DCM/ERF: Use fieldtrip functions to retrieve sensor location and
% orientation structure
hdr = read_ctf_res4(findres4file(S.Fdata));
grad = fieldtrip_ctf2grad(hdr);
D.channels.grad = grad;

% - coil locations (in this order - NZ:nazion , LE: left ear , RE: right ear)
CurrentDir = pwd;
cd(pre_data.folder);
hc_files = dir('*.hc');
if isempty(hc_files)
    warning(sprintf('Impossible to find head coil file\n'));
elseif length(hc_files) > 1
    hc_file = spm_select(1, '\.hc$', 'Select head coil file');
else
    hc_file = fullfile(pre_data.folder,hc_files.name);
end
clear hc_files
for coils=1:3
    fid = fopen(hc_file,'r');
    testlines = fgetl(fid);
    t=0;
    while t==0
        if coils==1 & strmatch('measured nasion coil position relative to head (cm):',testlines);
            t=1;
            for i = 1:3 % Nazion coordinates
                UsedLine    = fgetl(fid);
                UsedLine    = fliplr(deblank(fliplr(UsedLine)));
                [A,COUNT,ERRMSG,NEXTINDEX] = sscanf(UsedLine,'%c = %f');
                if ~isempty(ERRMSG) | (COUNT ~= 2)
                    warning(sprintf('Unable to read head coil file\n'));
                else
                    NZ(i) = A(2);
                end
            end
        end
        if coils==2 & strmatch('measured left ear coil position relative to head (cm):',testlines);
            t=1;
            for i = 1:3 % Nazion coordinates
                UsedLine    = fgetl(fid);
                UsedLine    = fliplr(deblank(fliplr(UsedLine)));
                [A,COUNT,ERRMSG,NEXTINDEX] = sscanf(UsedLine,'%c = %f');
                if ~isempty(ERRMSG) | (COUNT ~= 2)
                    warning(sprintf('Unable to read head coil file\n'));
                else
                    LE(i) = A(2);
                end
            end
        end
        if coils==3 & strmatch('measured right ear coil position relative to head (cm):',testlines);
            t=1;
            for i = 1:3 % Nazion coordinates
                UsedLine    = fgetl(fid);
                UsedLine    = fliplr(deblank(fliplr(UsedLine)));
                [A,COUNT,ERRMSG,NEXTINDEX] = sscanf(UsedLine,'%c = %f');
                if ~isempty(ERRMSG) | (COUNT ~= 2)
                    warning(sprintf('Unable to read head coil file\n'));
                else
                    RE(i) = A(2);
                end
            end
        end
        testlines = fgetl(fid);
    end
    fclose(fid);
end


CoiLoc = 10*[NZ ; LE ; RE]; % convertion from cm to mm
cd(CurrentDir);

fic_sensloc   = fullfile(D.path,[nam '_fidloc_meeg.mat']);
save(fic_sensloc,'CoiLoc');
clear hc_file CoiLoc UnusedLines UsedLine A COUNT ERRMSG NEXTINDEX
% -------

D.modality = 'MEG';
D.units = 'femto T';

if spm_matlab_version_chk('7') >= 0
    save(fullfile(D.path, D.fname), '-V6', 'D');
else
	save(fullfile(D.path, D.fname), 'D');
end

% find file name if truncated or with uppercase extension
% added by Arnaud Delorme June 15, 2004
% -------------------------------------------------------
function res4name = findres4file( folder )

res4name = dir([ folder filesep '*.res4' ]);
if isempty(res4name)
    res4name = dir([ folder filesep '*.RES4' ]);
end

if isempty(res4name)
    error('No file with extension .res4 or .RES4 in selected folder');
else
    res4name = [ folder filesep res4name.name ];
end;
return
