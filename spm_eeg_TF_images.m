function D = spm_eeg_TF_images(S)
% User interface for conversion of EEG-files to SPM's data structure
% FORMAT D = spm_eeg_converteeg2mat(S)
%
% struct S is optional and has the following (optional) fields:
%    fmt       - string that determines type of input file. Currently, this
%                string can be either 'CNT' or 'BDF'
%    Mname     - char matrix of input file name(s)
%    Fchannels - String containing name of channel template file
%_______________________________________________________________________
% 
% spm_eeg_converteeg2mat is a user interface to convert EEG-files from their
% native format to SPM's data format. This function assembles some
% necessary information before branching to the format-specific conversion
% routines.
% The user has to specify, by either using struct S or the GUI, a 'channel
% template file' that contains information about the (approximate) spatial
% positions of the channels.
% 
% Output: The converted data are written to files. The header
% structs, but not the data, are returned in D as a cell vector of structs.
%_______________________________________________________________________
%
% Additional formats can be added by (i) extending the code below in a
% straightforward fashion, (ii) providing a new channel template file and
% (iii) adding the actual conversion routine to the SPM-code.
%_______________________________________________________________________
% %W% Stefan Kiebel %E% 

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','TF',0);
try
	D = S.D;
catch
	D = spm_get(1, '.mat', 'Select EEG mat file');
end
P = spm_str_manip(D, 'H');

try
	D = spm_eeg_ldata(D);
catch    
	error(sprintf('Trouble reading file %s', D));
end

if isfield(D, 'Nfrequencies');
	try
		fmt = S.fmt;
	catch
		spm_input('average over ...', 1, 'd')
		Ctype = {
			'electrodes',...
				'frequency'};
		str   = 'Average over which dimension';
		Sel   = spm_input(str, 2, 'm', Ctype);
		fmt = Ctype{Sel};
	end
	
	switch fmt
		case {'electrodes'}
			try
				elecs = S.thresholds.elecs;
			catch 
				str = 'electrodes[s]';
				Ypos = -1;
				
				while 1
					if Ypos == -1   
						[D.electrodes_of_interest, Ypos] = spm_input(str, '+1', 'r', [], [1 Inf]);
					else
						D.electrodes_of_interest = spm_input(str, Ypos, 'r', [], [1 Inf]);
					end
					t=1:D.Nchannels;
					tmp=[];
					for en=D.electrodes_of_interest;
						if isempty(find(t==en))
							tmp=[tmp,en];
						end
					end
					if isempty(tmp) break, end
				end
				%[P, F] = fileparts(spm_str_manip(fname(:), 'r'));
				%[m, sta] = mkdir(P, spm_str_manip(Fname(:), 'tr'));
				%cd(fullfile(P, F));
				
				for i = 1 : D.events.Ntypes
					Itrials = find(D.events.code == D.events.types(i) & ~D.events.reject);
					dname = sprintf('trialtype%d', D.events.types(i));
					[m, sta] = mkdir(dname);
					cd(dname);
					
					for l = Itrials
						% if single trial data make new directory for single trials,
						% otherwise just write images to trialtype directory
						if D.Nevents ~= D.events.Ntypes
							% single trial data
							dname = sprintf('trial%d.img', l);
							fname = dname;
							[m, sta] = mkdir(dname);
							cd(dname);
						else
							fname = 'average.img';
						end
						data=squeeze(mean(D.data(D.electrodes_of_interest,:,:,i),1));	
						V.fname = fname;
						V.dim = [D.Nfrequencies D.Nsamples  1 ];
						V.dt=[spm_type('double') 0]; %%%check later with john
						V.mat = eye(4);
						V.pinfo = [1 0 0]';
						
						spm_write_vol(V, data); % d is data
					end
					
				end
			end
		case {'frequency'}
			try
				freqs = S.freqs;
			catch 
				str = 'Frequency window';
				Ypos = -1;
				
				while 1
					if Ypos == -1
						Ypos = '+1';
					end
					[D.Frequency_window, Ypos] = spm_input(str, Ypos, 'r', [], 2);
					inds=find(D.tf.frequencies>=D.Frequency_window(1) & D.tf.frequencies<=D.Frequency_window(2))
					if ~isempty(inds) break, end
					str = 'No data in range';
				end
			end
			data=squeeze(mean(D.data(:,inds,:,:),2));
			D.fnamedat = ['F' num2str(D.Frequency_window(1)) '_' num2str(D.Frequency_window(2)) '_' D.fnamedat];
			fpd = fopen(fullfile(P, D.fnamedat), 'w');
			for i=1:D.Nevents;
				D.scale.values(:, i) = spm_eeg_write(fpd, squeeze(data(:,:,i)), 2, D.datatype);
			end
			D.scale.dim=[1 3];
			D=rmfield(D,'Nfrequencies');
			D=rmfield(D,'tf');
			D.fname = ['F' num2str(D.Frequency_window(1)) '_' num2str(D.Frequency_window(2)) '_' D.fname];
			D.data = [];
			cd(D.path)
			if str2num(version('-release'))>=14
				save(fullfile(P, D.fname), '-V6', 'D');
			else
				save(fullfile(P, D.fname), 'D');
			end
			S.Fname=D.fname;
			spm_eeg_convertmat2ana(S);
			pause
	end
else
	spm_eeg_convertmat2ana;
end