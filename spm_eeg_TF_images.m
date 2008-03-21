function D = spm_eeg_TF_images(S)
% User interface for conversion of EEG-files to SPM's data structure
% FORMAT D = spm_eeg_TF_images(S)
%
% struct S is optional and has the following (optional) fields:
%    fmt       - string that determines type of input file. Currently, this
%                string can be 'electrodes' or 'frequency'
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
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% James Kilner, Stefan Kiebel 
% $Id: spm_eeg_TF_images.m 1237 2008-03-21 14:54:07Z stefan $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','TF',0);
try
    D = S.D;
catch
    D = spm_select(1, '\.mat$', 'Select EEG mat file');
    
end
P = spm_str_manip(D, 'H');

try
    D = spm_eeg_load(D);
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
                electrodes_of_interest = S.thresholds.elecs;
            catch 
                str = 'electrodes[s]';
                Ypos = -1;
                
                while 1
                    if Ypos == -1   
                        [electrodes_of_interest, Ypos] = spm_input(str, '+1', 'r', [], [1 Inf]);
                    else
                        electrodes_of_interest = spm_input(str, Ypos, 'r', [], [1 Inf]);
                    end
                    
                    
                    t = 1:D.nchannels;
                    
                    tmp = [];
                    for en = electrodes_of_interest;
                        if isempty(find(t == en))
                            tmp=[tmp, en];
                        end
                    end
                    if isempty(tmp) break, end
                end
            end
            
            try
                Nregion = S.region_no;
            catch 
                str = 'region number';
                Ypos = -1;
                
                while 1
                    if Ypos == -1   
                        [Nregion, Ypos] = spm_input(str, '+1', 'r', [], [1 Inf]);
                    else
                        Nregion = spm_input(str, Ypos, 'r', [], [1 Inf]);
                    end
                    if ~isempty(Nregion) break, end
                    str = 'No data';
                end
                
            end
            
            cl = cellstr(D{k}.conditionlabels);

            for i = 1 : D.nconditions
                Itrials = intersect(pickconditions(D{k}, cl(i)), find(~D{k}.reject))';
                
                cd(D.path)
                dname = sprintf('%dROI_TF_trialtype_%s', Nregion, D.conditionlabels(i));
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

                    dat = file_array(fname,...
                        [D.nfrequencies D.nsamples], 'FLOAT64');
                    dat(:,:) = squeeze(mean(D(electrodes_of_interest, :, :, i), 1));
                    
                    N = nifti;
                    N.dat = dat;
                    N.mat = eye(4);
                    N.mat_intent = 'Aligned';
                    create(N);
                end
                
            end
            
        case {'frequency'}
            try
                Frequency_window = S.freqs;
                Ypos = -1;
                while 1
                    if Ypos == -1
                        Ypos = '+1';
                    end
                    
                    inds = find(D.tf.frequencies >= Frequency_window(1) & D.tf.frequencies <= Frequency_window(2));
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
                    
                    inds = find(D.tf.frequencies >= Frequency_window(1) & D.tf.frequencies <= Frequency_window(2));
                    if ~isempty(inds) break, end
                    str = 'No data in range';
                end
            end
            
            % generate new meeg object with new filenames
            fnamedat = ['F' num2str(Frequency_window(1)) '_' num2str(Frequency_window(2)) '_' D.fnamedat];
            Dnew = newdata(D, fnamedat, [D.nchannels D.nsamples D.ntrials], dtype(D));
            
            Dnew = putdata(Dnew, 1:Dnew.nchannels, 1:Dnew.nsamples, Dnew.ntrials,... 
                squeeze(mean(D(:, inds, :, :), 2)));

           
            D=rmfield(D,'Nfrequencies');
            D=rmfield(D,'tf');

            save(D);
            S.D = fullfile(D.path, D.fname);
            
            try
                n = S.n;
            catch
                S.n = spm_input('Output image dimension', '+1', 'n', '32', 1);
                n = S.n;
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
