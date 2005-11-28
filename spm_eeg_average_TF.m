
function D=spm_eeg_average_TF(S)
%%% function to average induced TF data if standard average does not work because of out of memeory issues

% James Kilner
% $Id$

try
	D = S.D;
catch
    D = spm_select(1, '.*\.mat$', 'Select EEG mat file');
end
P = spm_str_manip(D, 'H');
try
	D = spm_eeg_ldata(D);
catch    
	error(sprintf('Trouble reading file %s', D));
end
try
	c = S.c;
catch
	% if there is no S.c, assume that user wants default average within
	% trial type
	c = eye(D.events.Ntypes);
end

d=zeros(size(D.data,1),size(D.data,2),size(D.data,3),D.events.Ntypes);
fh=fopen(fullfile(D.path,D.fnamedat),'r');

for n=1:D.Nevents
	[m,i]=find(D.events.types==D.events.code(n));
	data=fread(fh,[1,size(D.data,1)*size(D.data,2)*size(D.data,3)],'short');
	data=reshape(data,size(D.data,1),size(D.data,2),size(D.data,3),1);
	d(:,:,:,i)=d(:,:,:,i)+data;
end

D.fnamedat = ['m' D.fnamedat];

fpd = fopen(fullfile(P, D.fnamedat), 'w');
D.scale = zeros(D.Nchannels, 1, 1, D.events.Ntypes);
for n=1:D.events.Ntypes
	dat=squeeze(d(:,:,:,n)./length(find(D.events.code==D.events.types(n))));
	D.scale(:, 1, 1, n) = max(max(squeeze(abs(dat)), [], 3), [], 2)./32767;
	dat = int16(dat./repmat(D.scale.values(:, n), [1, D.Nfrequencies, D.Nsamples]));
	fwrite(fpd, dat, 'int16');
end
fclose (fh)	;

fclose(fpd);
D.Nevents = size(c, 2);

% labeling of resulting contrasts, take care to keep numbers of old trial
% types
% check this again: can be problematic, when user mixes within-trialtype
% and over-trial type contrasts
D.events.code = size(1, size(c, 2));
for i = 1:size(c, 2)
	if sum(c(:, i)) == 1 & sum(~c(:, i)) == size(c, 1)-1
		D.events.code(i) = find(c(:, i));
	else
		D.events.code(i) = i;
	end
end

D.events.time = [];
D.events.types = D.events.code;
D.events.Ntypes = length(D.events.types);
D.data = [];
D.events.reject = zeros(1, D.Nevents);
D.events.blinks = zeros(1, D.Nevents);
D.fname = ['m' D.fname];
if str2num(version('-release'))>=14
	save(fullfile(P, D.fname), '-V6', 'D');
else
	save(fullfile(P, D.fname), 'D');
end

