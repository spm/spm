%%%% function to read in CTF data to Matlab

pre_data=ctf_read;
PP1=squeeze(pre_data.data(:,1));
PP2=squeeze(pre_data.data(:,2));

inds=find(diff(PP1)>0);
D.events.code=PP1(inds+1)';
D.events.time=inds'+1;
inds=find(diff(PP2)<0);
D.events.code=[D.events.code,PP2(inds+1)'];
D.events.time=[D.events.time,inds'+1];

[X,I]=sort(D.events.time)
D.events.time=D.events.time(I);
D.events.code=D.events.code(I);


sens=pre_data.sensor.index.all_sens;
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

D.scale.dim = 1;
D.scale.values = zeros(D.Nchannels, 1);
fpd = fopen(fullfile(D.path, D.fnamedat), 'w');

fwrite(fpd, pre_data.data, 'float');
if str2num(version('-release'))>=14
    save(fullfile(D.path, D.fname), '-V6', 'D');
else
	save(fullfile(D.path, D.fname), 'D');
end


