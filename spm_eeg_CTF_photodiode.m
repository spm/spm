function D = spm_eeg_CTF_photodiode(S)
%%%% function to read in CTF data to Matlab
Mname = spm_select(inf, 'dir');
S.Fdata = deblank(Mname(1, 1:end-1));
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


I = STRMATCH('UPPT0',pre_data.sensor.label)
if isempty(I)
    error(sprintf('No parallel port channel was found in the CTF file'))
end
PP1=squeeze(pre_data.data(:,I(1)));
PP2=squeeze(pre_data.data(:,I(2)));

inds=find(diff(PP1)>0);
D.events.time=[];
D.events.code=[];
if ~isempty(inds)
    if length(PP1)<inds(end)+2
        inds(end)='';
    end
    D.events.code=PP1(inds+2)'; %changed to +2 from +1 to avoid errors when changing event code without passing by zero.
    D.events.time=inds'+1;
end

inds=find(diff(PP2)<0);
if ~isempty(inds)
    D.events.code=[D.events.code,PP2(inds+1)'+255];
    
    D.events.time=[D.events.time,inds'+1];
end
D.eventtypes=unique(D.events.code);

[X,I]=sort(D.events.time);
D.events.time=D.events.time(I);
D.events.code=D.events.code(I);
D.Radc=pre_data.setup.sample_rate;
I = STRMATCH('UADC0',pre_data.sensor.label)
if isempty(I)
    error(sprintf('No ADC channel was found in the CTF file'))
end
photodiode=pre_data.data(:,I);
[B,A] = BUTTER(3,40/(D.Radc/2));
photodiode=FILTFILT(B, A, photodiode);
clear pre_data;
xax=0:1:50;
xax=xax*1000/D.Radc;

for ty=1:length(D.eventtypes)
    times=D.events.time(D.events.code==D.eventtypes(ty));
    figure
    for n=1:length(times);
        plot(xax,photodiode(times(n):times(n)+50,:))
        hold on
        xlabel('Time (ms)')
    end
    title(['Event type ' num2str(D.eventtypes(ty))])
end
