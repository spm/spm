function Do = spm_eeg_grandmean(S)
% average over multiple data sets
% FORMAT Do = spm_eeg_grandmean(S)
% 
% S		    - struct (optional)
% (optional) fields of S:
% P			- filenames (char matrix) of EEG mat-file containing epoched
%             data  
% Pout      - filename (with or without path) of output file
% 
% Output:
% Do		- EEG data struct, result files are saved in the same
%                 directory as first input file.
%_______________________________________________________________________
% 
% spm_eeg_grandmean averages data over multiple files. The data must have
% the same trialtype numbering and sampling rate. This function can be used for
% grand mean averaging, i.e. computing the average over multiple subjects.
% Missing event types and bad channels are taken into account properly. 
% The output is written to a user-specified new file. The default name is 
% the same name as the first selected input file, but prefixed with a 'g'.
% The output file is written to the current working directory.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id: spm_eeg_grandmean.m 710 2006-12-21 14:59:04Z stefan $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG grandmean setup', 0);

try
    P = S.P;
catch
    P = spm_select(inf, '\.mat$', 'Select EEG mat files');
end

try
    S.Pout;
catch
    [filename, pathname] = uiputfile('*.mat', 'Select output file please');
    S.Pout = fullfile(pathname, filename);
end

D = cell(1, size(P, 1));
for i = 1:size(P, 1)
    try
        D{i} = spm_eeg_ldata(deblank(P(i, :)));
    catch    
        error(sprintf('Trouble reading files %s', deblank(P(i, :))));
    end
end

Nfiles = length(D);

% check dimensionality of data
dim = []; s = [];
for i = 1:Nfiles
    dim(i,:) = [length(D{i}.channels.eeg) size(D{i}.data, [2:length(size(D{i}.data))])];
    s(i) = D{i}.Radc;
    
    try
        D{i}.channels.Bad;
    catch
        D{i}.channels.Bad = [];
    end

end

if ~all(dim(:,1)==dim(1,1))
    ind = find(dim(1,1)~=dim(:,1));
    error(sprintf('Data don''t have the same number of channels.\nThere is a difference between files %s and %s.', D{1}.fname, D{ind(1)}.fname));
end

if ~all(dim(:,2)==dim(1,2))
    ind = find(dim(1,2)~=dim(:,2));
    error(sprintf('Data don''t have the same number of time points.\nThere is a difference between files %s and %s.', D{1}.fname, D{ind(1)}.fname));
end

if ~all(s == s(1))
    ind = find(s(1)~=s);
    error(sprintf('Data don''t have the same sampling rates.\nThere is a difference between files %s and %s.', D{1}.fname, D{ind(1)}.fname));
end


% output
Do = D{1};

try
    [f1, f2, f3] = fileparts(S.Pout);
    Do.path = f1;
    Do.fnamedat = [f2 '.dat'];    
catch
    Do.fnamedat = ['g' D{1}.fnamedat];
end

fpd = fopen(fullfile(Do.path, Do.fnamedat), 'w');

spm('Pointer', 'Watch'); drawnow;

if isfield(D{1}, 'Nfrequencies');
	Do.scale = zeros(length(D{1}.channels.eeg), 1, 1, D{1}.events.Ntypes);
	
	for i = 1:D{1}.events.Ntypes
		d = zeros(length(D{1}.channels.eeg), D{1}.Nfrequencies, D{1}.Nsamples);
		for k = 1:length(D)
			d = d+squeeze(D{k}.data(D{k}.channels.eeg,:,:,i));
		end
		d = d./length(D);
		
		Do.scale(:, 1, 1, i) = max(max(abs(d), [], 3), [], 2)./32767;
		d = int16(d./repmat(Do.scale(:, 1, 1, i), [1, Do.Nfrequencies, Do.Nsamples]));
		fwrite(fpd, d, 'int16');
        w=[];
	end
	
else

	Do.scale = zeros(length(D{1}.channels.eeg), 1, D{1}.Nevents);
	
    % how many different trial types and bad channels
    types = [];
    for i = 1:Nfiles
        types = unique([types D{i}.events.types]);       
    end

    Ntypes = length(types);
    
    % for determining bad channels of the grandmean
    w = zeros(length(D{1}.channels.eeg), Ntypes);
    
    % how many repetitons per trial type
    repl = zeros(1, Ntypes);
    for i = 1:Nfiles
        for j = 1:D{i}.events.Ntypes
            ind = find(D{i}.events.types(j) == types);
            repl(ind) =  repl(ind) + D{i}.events.repl(j);
        end
    end

    spm_progress_bar('Init', Ntypes, 'Events averaged'); drawnow;
    if Ntypes > 100, Ibar = floor(linspace(1, Ntypes, 100));
    else, Ibar = [1:Ntypes]; end

    for i = 1:Ntypes
        d = zeros(length(D{1}.channels.eeg), D{1}.Nsamples);

        for j = 1:length(D{1}.channels.eeg)

            for k = 1:Nfiles
                if ~ismember(j, D{k}.channels.Bad)
                    if ismember(types(i), D{k}.events.types)
                        d(j, :) = d(j, :) + D{k}.data(j, :, find(D{k}.events.types == types(i)));
                        w(j, i) = w(j, i) + 1;
                    end
                end
            end

            if w(j, i) > 0
                d(j, :) = d(j, :)/w(j, i);
            end

        end
		
		Do.scale(:, 1, i) = spm_eeg_write(fpd, d, 2, Do.datatype);
		
        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end

	end
end
fclose(fpd);

Do.data = [];

Do.fname = [spm_str_manip(Do.fnamedat, 'r') '.mat'];


% Do.Nchannels = length(D{1}.channels.eeg);
% Do.channels.name = D{1}.channels.name(D{1}.channels.eeg);
% 
% % remove eog and ref channels
% 
% eog = [];
% if isfield(Do.channels, 'heog')
%     eog = [eog Do.channels.heog];
%     Do.channels.heog = 0;
% end
% 
% if isfield(Do.channels, 'veog')
%     eog = [eog Do.channels.veog];    
%     Do.channels.veog = 0;
% end
% 
% Do.channels.order(eog) = [];
% 
% if isfield(Do.channels, 'reference')
%     
%     Do.channels.reference = 0;
% end
D = Do;

D.channels.Bad = find(~any(w'));
D.events.types = types;
D.events.Ntypes = Ntypes;
D.events.code = types;
D.events.repl = repl;
D.Nevents = Ntypes;
D.events.reject = zeros(1, Ntypes);
D.events.blinks = zeros(1, Ntypes);


if spm_matlab_version_chk('7') >= 0
    save(fullfile(D.path, D.fname), '-V6', 'D');
else
    save(fullfile(D.path, D.fname), 'D');
end


spm('Pointer', 'Arrow');
