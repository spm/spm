function D = spm_eeg_average(S);
% averages each channel over trials or trial types.
% FORMAT D = spm_eeg_average(S)
%
% S		    - optional input struct
% (optional) fields of S:
% D			- filename of EEG mat-file with epoched data
%
% Output:
% D			- EEG data struct (also written to files)
%_______________________________________________________________________
%
% spm_eeg_average averages single trial data within trial type. 
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_average.m 493 2006-04-04 20:02:47Z james $

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

D.fnamedat = ['m' D.fnamedat];

fpd = fopen(fullfile(P, D.fnamedat), 'w');

spm('Pointer', 'Watch'); drawnow;

if isfield(D, 'Nfrequencies');
    D.scale = zeros(D.Nchannels, 1, 1, D.events.Ntypes);

    for i = 1:D.events.Ntypes
        d = mean(D.data(:,:,:, find((D.events.code == D.events.types(i)) & ~D.events.reject)), 4);

        D.scale(:, 1, 1, i) = max(max(abs(d), [], 3), [], 2)./32767;
        d = int16(d./repmat(D.scale(:, 1, 1, i), [1, D.Nfrequencies, D.Nsamples]));
        fwrite(fpd, d, 'int16');
    end
    Ntypes = D.events.Ntypes;
else

    if isfield(D, 'weights');
        d = zeros(D.Nchannels, D.Nsamples);
        D.scale = zeros(D.Nchannels, 1, D.events.Ntypes);
        Ntypes = D.events.Ntypes;
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
					for nl=(find(D.events.code==D.events.types(i)));
						tempwf=[tempwf,D.weights(j,(nl-1)*D.Nsamples+1:nl*D.Nsamples)];
					end
                    data=data';
                    tempwf=reshape(tempwf,D.Nsamples,length(find(D.events.code==D.events.types(i))));
                    
                    for t=1:size(data,1)
                        B(t)=sum(tempwf(t,:).*data(t,:))/sum(tempwf(t,:));
                    end
                    d(j, :) =B;
				else
					d(j,:)=zeros(1,D.Nsamples);
				end
			end
			D.scale(:, 1, i) = spm_eeg_write(fpd, d, 2, D.datatype);
			
		end           
	else
		
        Ntypes = D.events.Ntypes;
		D.scale = zeros(D.Nchannels, 1, Ntypes);
		
		spm_progress_bar('Init', Ntypes, 'Averages done'); drawnow;
		if Ntypes > 100, Ibar = floor(linspace(1, Ntypes, 100));
		else, Ibar = [1:Ntypes]; end
		
		for i = 1:Ntypes
			
            w = (D.events.code == D.events.types(i) & ~D.events.reject)';
			
			ni(i) = length(find(w));
			
			d = zeros(D.Nchannels, D.Nsamples);
			
			if ni(i) == 0
				warning('%s: No trials for trial type %d', D.fname, D.events.types(i)); 
            else
                w = w./sum(w); % vector of trial-wise weights
				for j = 1:D.Nchannels
					d(j, :) = w'*squeeze(D.data(j, :, :))';
				end
            end
			
			D.scale(:, 1, i) = spm_eeg_write(fpd, d, 2, D.datatype);
			
			if ismember(i, Ibar)
				spm_progress_bar('Set', i);
				drawnow;
			end
			
		end
	end
	
	spm_progress_bar('Clear');
end
fclose(fpd);

D.Nevents = Ntypes;

D.events.code = D.events.types;
D.events.time = [];

if ~isfield(D, 'Nfrequencies') & ~isfield(D, 'weights');

    D.events.repl = ni;
	disp(sprintf('%s: Number of replications per contrast:', D.fname))
	s = [];
	for i = 1:Ntypes
		s = [s sprintf('average %d: %d trials', D.events.types(i), D.events.repl(i))];
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
if isfield(D, 'weights');
    rmfield(D, 'weights');
end
D.fname = ['m' D.fname];

if str2num(version('-release'))>=14
	save(fullfile(P, D.fname), '-V6', 'D');
else
	save(fullfile(P, D.fname), 'D');
end

spm('Pointer', 'Arrow');
