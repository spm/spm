function D = spm_eeg_convert2images(S)
% User interface for conversion of EEG-files to SPM's data structure
% FORMAT D = spm_eeg_convert2images(S)
%
% struct S is optional and has the following (optional) fields:
%    fmt       - string that determines type of input file. Currently, this
%                string can be 'electrodes' or 'frequency'
%    Mname     - char matrix of input file name(s)
%    Fchannels - String containing name of channel template file
%_______________________________________________________________________
% 
% spm_eeg_convert2images is a user interface to convert EEG-files from their
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
% $Id: spm_eeg_convert2images.m 1278 2008-03-28 18:38:11Z stefan $

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

if strcmp(D.type, 'continuous')
    error('Data are continuous. Try epoched data.');
end

if strcmp(D.transformtype, 'TF');
    try
        images.fmt = S.images.fmt;
    catch
        spm_input('average over ...', 1, 'd')
        Ctype = {
            'electrodes',...
                'frequency'};
        str   = 'Average over which dimension';
        Sel   = spm_input(str, 2, 'm', Ctype);
        images.fmt = Ctype{Sel};
    end
    
    switch images.fmt
        case {'electrodes'}
            try
                images.electrodes_of_interest = S.images.elecs;
            catch 
                str = 'electrodes[s]';
                Ypos = -1;
                
                while 1
                    if Ypos == -1   
                        [images.electrodes_of_interest, Ypos] = spm_input(str, '+1', 'r', [], [1 Inf]);
                    else
                        images.electrodes_of_interest = spm_input(str, Ypos, 'r', [], [1 Inf]);
                    end
                    
                    if any(ismember(images.electrodes_of_interest, [1:D.nchannels]))
                        break
                    end
                end
            end
            
            try
                images.Nregion = S.images.region_no;
            catch 
                str = 'region number';
                Ypos = -1;
                
                while 1
                    if Ypos == -1   
                        [images.Nregion, Ypos] = spm_input(str, '+1', 'r', [], [1 Inf]);
                    else
                        images.Nregion = spm_input(str, Ypos, 'r', [], [1 Inf]);
                    end
                    if ~isempty(images.Nregion) break, end
                    str = 'No data';
                end
                
            end
            
            cl = unique(D.conditions);

            for i = 1 : D.nconditions
                Itrials = intersect(pickconditions(D, cl{i}), find(~D.reject))';
                
                cd(D.path)
                dname = sprintf('%dROI_TF_trialtype_%s', images.Nregion, cl{i});
                [m, sta] = mkdir(dname);
                cd(dname);
                
                for l = Itrials

					% if single trial data make new directory for single trials,
					% otherwise just write images to trialtype directory
                    if strcmp(D.type, 'single')
						% single trial data
						fname = sprintf('trial%d.img', l);
					else
						fname = 'average.img';
					end

                    dat = file_array(fname, [D.nfrequencies D.nsamples], 'FLOAT64');
                    dat(:, :) = squeeze(mean(D(images.electrodes_of_interest, :, :, l), 1));
                    
                    N = nifti;
                    N.dat = dat;
                    N.mat = [1 0 0  min(D.frequencies);
                        0 1000/D.fsample 0           time(D, 1, 'ms');
                        0 0 1           0;
                        0 0 0           1];
                    N.mat_intent = 'Aligned';
                    create(N);
                    
                end
                
            end
            
        case {'frequency'}
            try
                images.Frequency_window = S.images.freqs;
                Ypos = -1;
                while 1
                    if Ypos == -1
                        Ypos = '+1';
                    end
                    
                    inds = find(D.frequencies >= Frequency_window(1) & D.frequencies <= Frequency_window(2));
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
                    [images.Frequency_window, Ypos] = spm_input(str, Ypos, 'r', [], 2);
                    
                    inds = find(D.frequencies >= images.Frequency_window(1) & D.frequencies <= images.Frequency_window(2));
                    if ~isempty(inds) break, end
                    str = 'No data in range';
                end
            end
            
            % generate new meeg object with new filenames
            fnamedat = ['F' num2str(images.Frequency_window(1)) '_' num2str(images.Frequency_window(2)) '_' D.fnamedat];
            Dnew = clone(D, fnamedat, [D.nchannels D.nsamples D.ntrials]);
            
            Dnew(1:Dnew.nchannels, 1:Dnew.nsamples, 1:Dnew.ntrials) = ... 
                squeeze(mean(D(:, inds, :, :), 2));

           % fake time-series data
           Dnew = transformtype(Dnew, 'time');
           

            save(Dnew);
            S.Fname = fullfile(Dnew.path, Dnew.fname);
            
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
            spm_eeg_convertmat2nifti3D(S);
    end
else
    clear S;
    S.Fname = fullfile(D.path, D.fname);
    spm_eeg_convertmat2nifti3D(S);
end
