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
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% James Kilner 
% $Id: spm_eeg_TF_images.m 539 2006-05-19 17:59:30Z Darren $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','TF',0);
try
	D = S.D;
catch
	D = spm_select(1, '\.mat$', 'Select EEG mat file');
	
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
				D.electrodes_of_interest = S.thresholds.elecs;
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
			end
			try
				D.Nregion = S.region_no;
			catch 
				str = 'region number';
				Ypos = -1;
				
				while 1
					if Ypos == -1   
						[D.Nregion, Ypos] = spm_input(str, '+1', 'r', [], [1 Inf]);
					else
						D.Nregion = spm_input(str, Ypos, 'r', [], [1 Inf]);
					end
					if ~isempty(D.Nregion) break, end
					str = 'No data';
				end
				
			end
			
			for i = 1 : D.events.Ntypes
				Itrials = find(D.events.code == D.events.types(i) & ~D.events.reject);
				cd(D.path)
				dname = sprintf('%dROI_TF_trialtype%d', D.Nregion, D.events.types(i));
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
					V.dt=[spm_type('float64') 0]; %%%check later with john
					V.mat = eye(4);
					V.pinfo = [1 0 0]';
					
					spm_write_vol(V, data); % d is data
				end
				
			end
			
        case {'frequency'}
            try
                D.Frequency_window = S.freqs;
                Ypos = -1;
                while 1
                    if Ypos == -1
                        Ypos = '+1';
                    end
                    
                    inds=find(D.tf.frequencies>=D.Frequency_window(1) & D.tf.frequencies<=D.Frequency_window(2))
                    if ~isempty(inds) break, end
                    str = 'No data in range';
                end
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
                D.scale(:,1,i) = spm_eeg_write(fpd, squeeze(data(:,:,i)), 2, D.datatype);
            end
            fclose(fpd);
           
            D=rmfield(D,'Nfrequencies');
            D=rmfield(D,'tf');
            D.fname = ['F' num2str(D.Frequency_window(1)) '_' num2str(D.Frequency_window(2)) '_' D.fname];
            D.data = [];
            cd(D.path)
            if spm_matlab_version_chk('7') >= 0,
                save(fullfile(P, D.fname), '-V6', 'D');
            else
                save(fullfile(P, D.fname), 'D');
            end
            S.Fname=D.fname;
            
            try
                n = S.n;
            catch
                S.n = spm_input('Output image dimension', '+1', 'n', '32', 1);
                n=S.n;
            end
            
            if length(n) > 1
                error('n must be scalar');
            end
            
            try
                interpolate_bad = S.interpolate_bad;
            catch
                S.interpolate_bad = spm_input('Interpolate bad channels or mask out?',...
                    '+1', 'b', 'Interpolate|Mask out', [1,0]);
            end
            spm_eeg_convertmat2ana(S);
    end
else
    clear S;
    S.Fname = fullfile(D.path, D.fname);
    spm_eeg_convertmat2ana(S);
end
