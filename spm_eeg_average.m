function D = spm_eeg_average(S);
% averages each channel over trials or trial types.
% FORMAT D = spm_eeg_average(S)
%
% S		    - optional input struct
% (optional) fields of S:
% D			- filename of EEG mat-file with epoched data
% c         - contrast matrix, each column computes a contrast of the data
%
% Output:
% D			- EEG data struct (also written to files)
%_______________________________________________________________________
%
% spm_eeg_average averages data, either over single trials or averaged
% data. By default, if applied to single trial data, the function will
% average within trial type. Additionally, one can supply a contrast matrix
% c that can be used to average over trial types or take differences
% between trial types.
%
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_average.m 112 2005-05-04 18:20:52Z john $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG averaging setup',0);

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

if size(c, 1) ~= D.events.Ntypes
	error('length of contrast must be number of trial types');
end

if sign(c(:)) ~= c(:)
	error('contrast can only take 1, -1 or 0');
end


D.fnamedat = ['m' D.fnamedat];

fpd = fopen(fullfile(P, D.fnamedat), 'w');

spm('Pointer', 'Watch'); drawnow;


if isfield(D, 'Nfrequencies');
	D.scale.dim = [1 4];
	D.scale.values = zeros(D.Nchannels, D.events.Ntypes);
	
	for i = 1:D.events.Ntypes
		d = mean(D.data(:,:,:, find((D.events.code == D.events.types(i)) & ~D.events.reject)), 4);
		
		D.scale.values(:, i) = max(max(abs(d), [], 3), [], 2)./32767;
		d = int16(d./repmat(D.scale.values(:, i), [1, D.Nfrequencies, D.Nsamples]));
		fwrite(fpd, d, 'int16');
	end
	
else
	
	if isfield(D, 'weights');	
		d = zeros(D.Nchannels, D.Nsamples);
		D.scale.dim = [1 3];
		D.scale.values = zeros(D.Nchannels, D.events.Ntypes);
		for i = 1:D.events.Ntypes
			
			for j = 1:D.Nchannels
				tempwf=[];
				ti=0;
				ts=0;
				while ts==0
					ti=ti+1;
					ts=(j==D.channels.thresholded{ti});
					
				end
				
				if isempty(ts)
					data=squeeze(D.data(j,:,find(D.events.code==D.events.types(i))))';
					ndata=reshape(data',size(data,1)*size(data,2),1);
					Xs=sparse(repmat(speye(size(data,2)),[size(data,1),1]));
					Wis=speye(length(ndata));
					for nl=(find(D.events.code==D.events.types(i)));
						
						tempwf=[tempwf,D.weights(j,(nl-1)*D.Nsamples+1:nl*D.Nsamples)];
					end
					Wis=spdiags(tempwf',0,Wis);
					
					d(j, :) =(((Xs'*Wis*Xs)^-1)*Xs'*Wis*ndata)';
				else
					d(j,:)=zeros(1,D.Nsamples);
				end
			end
			D.scale.values(:, i) = spm_eeg_write(fpd, d, 2, D.datatype);
			
		end           
		
		for i = 1:size(c, 2)
			
			c_tmp = zeros(D.Nevents, 1);
			for j = 1:size(c, 1)
				if c(j, i)
					ind = (D.events.code == D.events.types(j) & ~D.events.reject)';
					c_tmp = c_tmp + c(j, i)*ind;
				end
			end
			
			% weight contrast vectors
			if ~isempty(find(c_tmp == 1))
				w1 = sum(c_tmp == 1);
				c_tmp(c_tmp == 1) = 1/w1*c_tmp(c_tmp == 1);
			end
			if ~isempty(find(c_tmp == -1))
				wm1 = sum(c_tmp == -1);
				c_tmp(c_tmp == -1) = 1/wm1*c_tmp(c_tmp == -1);
			end
			
			% identify the case when there is only one observation (single
			% trial)
			ni(i) = sum(c_tmp ~= 0);
		end
		
		
	else
		
		D.scale.dim = [1 3];
		D.scale.values = zeros(D.Nchannels, size(c, 2));
		
		spm_progress_bar('Init', size(c, 2), 'Averages done'); drawnow;
		if size(c, 2) > 100, Ibar = floor(linspace(1, size(c, 2),100));
		else, Ibar = [1:size(c, 2)]; end
		
		for i = 1:size(c, 2)
			
			c_tmp = zeros(D.Nevents, 1);
			for j = 1:size(c, 1)
				if c(j, i)
					ind = (D.events.code == D.events.types(j) & ~D.events.reject)';
					c_tmp = c_tmp + c(j, i)*ind;
				end
			end
			
			% weight contrast vectors
			if ~isempty(find(c_tmp == 1))
				w1 = sum(c_tmp == 1);
				c_tmp(c_tmp == 1) = 1/w1*c_tmp(c_tmp == 1);
			end
			if ~isempty(find(c_tmp == -1))
				wm1 = sum(c_tmp == -1);
				c_tmp(c_tmp == -1) = 1/wm1*c_tmp(c_tmp == -1);
			end
			
			% identify the case when there is only one observation (single
			% trial)
			ni(i) = sum(c_tmp ~= 0);
			
			d = zeros(D.Nchannels, D.Nsamples);
			
			if ni(i) == 0
				warning('%s: No trials for contrast %d', D.fname, i); 
			else
				for j = 1:D.Nchannels
					d(j, :) = c_tmp'*squeeze(D.data(j, :, :))';           
					
				end
			end           
			
			D.scale.values(:, i) = spm_eeg_write(fpd, d, 2, D.datatype);
			
			if ismember(i, Ibar)
				spm_progress_bar('Set', i);
				drawnow;
			end
			
		end
	end
	
	spm_progress_bar('Clear');
end
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

if ~isfield(D, 'Nfrequencies');
	D.events.repl = ni;
	disp(sprintf('%s: Number of replications per contrast:', D.fname))
	s = [];
	for i = 1:size(c, 2)
		s = [s sprintf('Contrast %d: %d', D.events.types(i), D.events.repl(i))];
		if i < D.events.Ntypes
			s = [s sprintf(', ')];
		else
			s = [s '\n'];
		end
	end 
	disp(sprintf(s))
end

D.data = [];
D.events.reject = zeros(1, D.Nevents);
D.events.blinks = zeros(1, D.Nevents);

D.fname = ['m' D.fname];

if str2num(version('-release'))>=14
	save(fullfile(P, D.fname), '-V6', 'D');
else
	save(fullfile(P, D.fname), 'D');
end

spm('Pointer', 'Arrow');
