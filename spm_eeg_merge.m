function Dout = spm_eeg_merge(S)
% Concatenate epoched single trial files.
% FORMAT D = spm_eeg_merge(S)
%
% S           - input structure (optional)
%  fields of S:
%   S.D       - character array containing filename of M/EEG mat-files
%               or cell array of D's
%   S.recode  - this field specifies how the condition labels will be
%               translated from the original files to the merged file.
%               Several options are possible:
%                 'same'        - leave the condition labels unchanged
%                 'addfilename' - add the original file name to condition
%                                 label
%                 old way specification - (for backward compatibility)                     
%                       a cell array where each cell contains a condition
%                       label. The ordering of these labels must be such 
%                       that each row in the cell matrix specifies the 
%                       conditionlabels for one of the selected files.
%                 specification via recoding rules - for this S.recode
%                       should be a structure array where each element 
%                       specifies a rule using the following fields:
%                            file - can be a cell array of strings with 
%                                   file names, a vector of file indices 
%                                   or a string with regular expression 
%                                   matching the files to which the rule 
%                                   will apply.
%                            labelorg - can be a cell array of condition 
%                                   labels or a string with regular 
%                                   expression matching the condition 
%                                   labels to which this rule will apply.
%                            labelnew - new label for the merged file. It
%                                   can contain special tokens #file# and
%                                   #labelorg# that will be replaced by 
%                                   the original file name and original 
%                                   condition label respectively.
%                       The rule will be applied one after the other so 
%                       the last rule takes precedences. Trials not 
%                       matched by any of the rules will keep their 
%                       original labels.
%                       Example:
%                          S.recode(1).file     = '.*';
%                          S.recode(1).labelorg = '.*';
%                          S.recode(1).labelnew = '#labelorg# #file#';
%                       has the same effect as the 'addfilename' option.
%   S.prefix     - prefix for the output file (default - 'c')
%
% 
% Dout        - MEEG object (also written to disk)
%__________________________________________________________________________
%
% This function can be used to merge M/EEG files to one file. This is
% useful whenever the data are distributed over multiple files, but one
% wants to use all information in one file. For example, when displaying
% data (SPM displays data from only one file at a time), or merging
% information that has been measured in multiple sessions.
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging
%
% Stefan Kiebel, Vladimir Litvak, Doris Eckstein, Rik Henson
% $Id: spm_eeg_merge.m 6622 2015-12-03 11:54:13Z vladimir $

SVNrev = '$Rev: 6622 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG Merge'); spm('Pointer','Watch');

if ~isfield(S, 'prefix'),       S.prefix = 'c';           end
if ~isfield(S, 'recode'),       S.recode = 'same';        end

%-Load MEEG data
%--------------------------------------------------------------------------
D = S.D;

if ischar(D)
    F = cell(1,size(D,1));
    try
        for i = 1:size(D, 1)
            F{i} = spm_eeg_load(deblank(D(i, :)));
        end
        D = F;
    catch
        error('Trouble reading files');
    end
end

Nfiles = length(D);

if Nfiles < 2
    error('Need at least two files for merging');
end

%-Check input and determine number of new number of trial types
%--------------------------------------------------------------------------
Ntrials = [];
megsens = [];
eegsens = [];
fid     = [];
isTF    =  strncmpi(D{1}.transformtype,'TF',2); % TF and TFphase

for i = 1:Nfiles
    if ~isequal(D{i}.transformtype, D{1}.transformtype)
        error(['The datasets do not contain the same kind of data.\n'...
               'There is a difference between files\n\t%s\nand\n\t%s.'], ...
               D{1}.fname, D{i}.fname);
    end

    if D{1}.nchannels ~= D{i}.nchannels
        error(['Data don''t have the same number of channels.\n' ...
               'There is a difference between files\n\t%s\nand\n\t%s.'], ...
               D{1}.fname, D{i}.fname);
    end

    if D{1}.nsamples ~= D{i}.nsamples
        error(['Data don''t have the same number of time points.\n' ...
               'There is a difference between files\n\t%s\nand\n\t%s.'], ...
               D{1}.fname, D{i}.fname);
    end

    if D{1}.fsample ~= D{i}.fsample
        error(['Data don''t have the same sampling rate.\n' ...
               'There is a difference between files\n\t%s\nand\n\t%s.'], ...
               D{1}.fname, D{i}.fname);
    end

    if isTF &&  ~isequal(D{1}.frequencies, D{i}.frequencies)
        error(['Data don''t have the same frequencies.\n' ...
               'There is a difference between files\n\t%s\nand\n\t%s.'], ...
               D{1}.fname, D{i}.fname);
    end

    if ~isempty(D{i}.sensors('MEG'))
        megsens = spm_cat_struct(megsens, D{i}.sensors('MEG'));
    end
    
    if ~isempty(D{i}.sensors('EEG'))
        eegsens = spm_cat_struct(eegsens, D{i}.sensors('EEG'));
    end
    
    if ~isempty(megsens) || ~isempty(eegsens)
        fid = spm_cat_struct(fid, D{i}.fiducials);
    end
    
    Ntrials = [Ntrials D{i}.ntrials];
end

%-Prepare some useful lists
%--------------------------------------------------------------------------
F        = {};
Find     = [];
clb      = {};
for i = 1:Nfiles
    F{i} = fname(D{i});
    clb  = [clb D{i}.conditions];
    Find = [Find i*ones(1, D{i}.ntrials)];
end
uclb      = unique(clb);


%-Specify condition labels recoding
%--------------------------------------------------------------------------
if ~isfield(S, 'recode')
    S.recode = spm_input('What to do with condition labels?', 1, 'm',...
        'Leave as they are|Add file name|Specify rules for recoding|Specify recoding the old way', strvcat('same', 'addfilename', 'rules', 'old'));
end

if isequal(S.recode, 'old')
    S.recode = {};
    for i = 1:Nfiles
        for j = 1:nconditions(D{i})
            S.recode{i}{j} = spm_input(sprintf('Labels: %s', spm_file(D{i}.fname, 'basename')),...
                '+1', 's', D{i}.condlist{j});
        end
    end
elseif isequal(S.recode, 'rules')
    S.recode = [];
    stop     = 0;
    ind      = 1;   
    while ~stop
        spm_input(['Please define rule ' num2str(ind) ':'], 1, 'd');
        switch spm_input('To which files will this rule apply?', '+1', 'm',...
                'All the files|Specify indices|Select files|Wildcard expression (*,?)|Regular expression',...
                strvcat('all', 'indices', 'select', 'wildcard', 'regexp'))
            case 'all'
                S.recode(ind).file = '.*';
            case 'indices'
                S.recode(ind).file =  spm_input('Input file indices', '+1', 'n', num2str(ind), [1 Inf]);
            case 'select'
                [selection, ok]= listdlg('ListString', F, 'SelectionMode', 'multiple' ,'Name', 'Select files' , 'ListSize', [400 300]);
                if ok
                    S.recode(ind).file = F(selection);
                else
                    continue;
                end
            case 'wildcard'
                S.recode(ind).file = regexptranslate('wildcard' , spm_input('Input wildcard expresssion', '+1', 's',  '*'));
            case 'regexp'
                S.recode(ind).file = spm_input('Input regular expresssion', '+1', 's',  '.*');
        end

        switch spm_input('What conditions will be renamed?', '+1', 'm',...
                'All|Select|Specify by wildcard expression (*,?)|Specify by regular expression',...
                strvcat('all', 'select', 'wildcard', 'regexp'))
            case 'all'
                S.recode(ind).labelorg = '.*';
            case 'select'
                [selection, ok]= listdlg('ListString', uclb, 'SelectionMode', 'multiple' ,'Name', 'Select conditions' , 'ListSize', [400 300]);
                if ok
                    S.recode(ind).labelorg = uclb(selection);
                else
                    continue;
                end
            case 'wildcard'
                S.recode(ind).labelorg = regexptranslate('wildcard' , spm_input('Input wildcard expresssion', '+1', 's',  '*'));
            case 'regexp'
                S.recode(ind).labelorg = spm_input('Input regular expresssion', '+1', 's',  '.*');
        end

        S.recode(ind).labelnew = spm_input('Input the new name?', '+1', 's',  '');
     
        stop = spm_input('Define another rule?','+1','yes|stop', [0 1], 0);
        ind  = ind+1;
    end
end

%-Generate new meeg object with new filenames
%--------------------------------------------------------------------------
Dout = D{1};
[p, f, x] = fileparts(fnamedat(Dout));

if ~isTF
    Dout = clone(Dout, fullfile(pwd, [S.prefix f x]), [Dout.nchannels Dout.nsamples sum(Ntrials)]);
else
    Dout = clone(Dout, fullfile(pwd, [S.prefix f x]), [Dout.nchannels Dout.nfrequencies Dout.nsamples sum(Ntrials)]);
end


%-Perform condition labels recoding
%--------------------------------------------------------------------------
if isequal(S.recode, 'same')
    Dout = conditions(Dout, ':', clb);
elseif isequal(S.recode, 'addfilename')
    for i = 1:numel(clb)
        clb{i} = [clb{i} ' ' spm_file(F{Find(i)}, 'basename')];
    end
    Dout = conditions(Dout, ':', clb);
elseif iscell(S.recode)
    for i = 1:Nfiles
        ind = find(Find == i);        
               
        for j = 1:D{i}.nconditions
            clb(ind(strmatch(D{i}.condlist{j}, clb(ind), 'exact'))) = S.recode{i}(j);
        end
    end
    Dout = conditions(Dout, ':', clb);
elseif isstruct(S.recode)
    clbnew = clb;
    
    for i = 1:numel(S.recode)
        if isnumeric(S.recode(i).file)
            ind = S.recode(i).file;
        elseif iscell(S.recode(i).file)
            ind = spm_match_str(F, S.recode(i).file);
        elseif ischar(S.recode(i).file)
            ind = find(~cellfun('isempty', regexp(F, S.recode(i).file)));
        else
            error('Invalid file specification in recoding rule.');
        end
        
        ind = find(ismember(Find, ind));
        
        if iscell(S.recode(i).labelorg)
            ind = ind(ismember(clb(ind), S.recode(i).labelorg));
        elseif ischar(S.recode(i).labelorg)
            ind = ind(~cellfun('isempty', regexp(clb(ind), S.recode(i).labelorg)));
        else
            error('Invalid original condition label specification in recoding rule.');
        end
        
        for j = 1:length(ind)
            labelnew = S.recode(i).labelnew;
            labelnew = strrep(labelnew, '#file#', spm_file(F{Find(ind(j))}, 'basename'));
            labelnew = strrep(labelnew, '#labelorg#', clb(ind(j)));
            
            clbnew{ind(j)} = labelnew;
        end
    end
    
    Dout = conditions(Dout, ':', clbnew);
end
            
%-Average sensor locations
%--------------------------------------------------------------------------
if ~isempty(megsens)
    spm_figure('GetWin','Graphics');clf;
    if ~isempty(eegsens)
        h = subplot(2, 1, 1);
        aeegsens = ft_average_sens(eegsens, 'weights', Ntrials, 'feedback', h);
        Dout = sensors(Dout, 'EEG', aeegsens);
        
        h = subplot(2, 1, 2);
    else
        h = axes;
    end
    
    [amegsens,afid] = ft_average_sens(megsens, 'fiducials', fid, 'weights', Ntrials, 'feedback', h);
    Dout = sensors(Dout, 'MEG', amegsens);
    Dout = fiducials(Dout, afid);
elseif ~isempty(eegsens)
    spm_figure('GetWin','Graphics');clf;
    h = axes;
    [aeegsens,afid] = ft_average_sens(eegsens, 'fiducials', fid, 'weights', Ntrials, 'feedback', h);
    Dout = sensors(Dout, 'EEG', aeegsens);
    Dout = fiducials(Dout, afid);
end



%-Write files
%--------------------------------------------------------------------------
spm_progress_bar('Init', Nfiles, 'Files merged');
if Nfiles > 100, Ibar = floor(linspace(1, Nfiles,100));
else Ibar = [1:Nfiles]; end

k = 0;

for i = 1:Nfiles

    ind = union(Dout.badchannels, D{i}.badchannels);
    if ~isempty(ind)
        Dout = badchannels(Dout, ind, 1);
    end

    % write trial-wise to avoid memory mapping error
    for j = 1:D{i}.ntrials
        k = k + 1;
        if ~isTF
            Dout(1:Dout.nchannels, 1:Dout.nsamples, k) =  D{i}(:,:,j);
        else
            Dout(1:Dout.nchannels, 1:Dout.nfrequencies, 1:Dout.nsamples, k) =  D{i}(:,:,:,j);
        end
        Dout = badtrials(Dout, k, badtrials(D{i}, j));
    end
    
    % Propagate some useful information from the original files to the
    % merged file
    Dout = repl(Dout, find(Find == i), D{i}.repl);
    Dout = trialonset(Dout, find(Find == i), D{i}.trialonset);
    Dout = trialtag(Dout, find(Find == i), D{i}.trialtag);
    Dout = events(Dout, find(Find == i), D{i}.events);
    
    if ismember(i, Ibar), spm_progress_bar('Set', i); end

end

%-Save new M/EEG data
%--------------------------------------------------------------------------
Dout = Dout.history('spm_eeg_merge', S, 'reset');
save(Dout);

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
spm('FigName','M/EEG merge: done'); spm('Pointer','Arrow');
