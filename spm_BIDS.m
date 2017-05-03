function varargout = spm_BIDS(varargin)
% Parse directory structure formated according to the BIDS standard
% FORMAT BIDS = spm_BIDS(root)
% root   - directory formated according to BIDS [Default: pwd]
% BIDS   - structure containing the BIDS file layout
%
% FORMAT result = spm_BIDS(BIDS,query,...)
% BIDS   - BIDS directory name or structure containing the BIDS file layout
% query  - type of query
% result - query's result
%__________________________________________________________________________
%
% BIDS (Brain Imaging Data Structure): http://bids.neuroimaging.io/
%   The brain imaging data structure, a format for organizing and
%   describing outputs of neuroimaging experiments.
%   K. J. Gorgolewski et al, Scientific Data, 2016.
%__________________________________________________________________________
% Copyright (C) 2016-2017 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_BIDS.m 7073 2017-05-03 14:36:33Z guillaume $


%-Validate input arguments
%==========================================================================
if ~nargin
    root = pwd;
elseif nargin == 1
    if ischar(varargin{1})
        root = spm_select('CPath',varargin{1});
    else
        varargout = varargin(1);
        return;
    end
else
    BIDS = spm_BIDS(varargin{1});
    varargout{1} = BIDS_query(BIDS,varargin{2:end});
    return;
end

%-BIDS structure
%==========================================================================

BIDS = struct(...
    'dir',root, ...               % BIDS directory
    'description',struct([]), ... % content of dataset_description.json
    'sessions',{{}},...           % cellstr of sessions
    'scans',struct([]),...        % content of sub-<participant_label>_scans.tsv (should go within subjects)
    'sess',struct([]),...         % content of sub-participants_label>_sessions.tsv (should go within subjects)
    'participants',struct([]),... % content of participants.tsv
    'subjects',struct([]));       % structure array of subjects

%-Validation of BIDS root directory
%==========================================================================
if isempty(BIDS.dir)
    error('A BIDS directory has to be specified.');
elseif ~exist(BIDS.dir,'dir')
    error('BIDS directory does not exist.');
elseif ~exist(fullfile(BIDS.dir,'dataset_description.json'),'file')
    error('BIDS directory not valid: missing dataset_description.json.');
end

%-Dataset description
%==========================================================================
try
    BIDS.description = spm_jsonread(fullfile(BIDS.dir,'dataset_description.json'));
    if ~isfield(BIDS.description,'BIDSVersion') || ~isfield(BIDS.description,'Name')
        error('BIDS dataset description not valid.');
    end
catch
    error('BIDS dataset description could not be read.');
end
% See also optional README and CHANGES files

%-Optional directories
%==========================================================================
% [code/]
% [derivatives/]
% [stimuli/]
% [sourcedata/]

%-Scans key file
%==========================================================================

% sub-<participant_label>/[ses-<session_label>/]
%     sub-<participant_label>_scans.tsv

%-Participant key file
%==========================================================================
p = spm_select('FPList',BIDS.dir,'^participants\.tsv$');
if ~isempty(p)
    BIDS.participants = spm_load(p);
end
p = spm_select('FPList',BIDS.dir,'^participants\.json$');
if ~isempty(p)
    spm_jsonread(p); % need to be stored somewhere
end

p = spm_select('FPList',fullfile(BIDS.dir,'phenotype'),'.*\.tsv$');
p = spm_select('FPList',fullfile(BIDS.dir,'phenotype'),'.*\.json$');
% <measurement_tool_name>

%-Sessions file
%==========================================================================

% sub-<participant_label>/[ses-<session_label>/]
%      sub-<participant_label>[_ses-<session_label>]_sessions.tsv

%-Tasks: JSON files are accessed through metadata
%==========================================================================
%t = spm_select('FPList',BIDS.dir,...
%    '^task-.*_(beh|bold|events|channels|physio|stim|meg)\.(json|tsv)$');

%-Subjects
%==========================================================================
sub = cellstr(spm_select('List',BIDS.dir,'dir','^sub-.*$'));
if isequal(sub,{''})
    error('No subjects found in BIDS directory.');
end

for su=1:numel(sub)
    sess = cellstr(spm_select('List',fullfile(BIDS.dir,sub{su}),'dir','^ses-.*$'));    
    for se=1:numel(sess)
        if isempty(BIDS.subjects)
            BIDS.subjects = parse_subject(BIDS.dir, sub{su}, sess{se});
        else
            BIDS.subjects(end+1) = parse_subject(BIDS.dir, sub{su}, sess{se});
        end
    end
end

varargout = { BIDS };


%==========================================================================
%-Parse a subject's directory
%==========================================================================
function subject = parse_subject(p, subjname, sesname)

subject.name = subjname; % subject name ('sub-<participant_label>')
subject.path = fullfile(p,subjname,sesname); % full path to subject directory
subject.session = sesname;  % session name ('' or 'ses-*')
subject.anat = struct([]);  % anatomy imaging data
subject.func = struct([]);  % task imaging data
subject.fmap = struct([]);  % fieldmap data
subject.beh = struct([]);   % behavioral experiment data
subject.dwi = struct([]);   % diffusion imaging data
subject.meg = struct([]);   % MEG data
subject.pet = struct([]);   % PET imaging data


%--------------------------------------------------------------------------
%-Anatomy imaging data
%--------------------------------------------------------------------------
pth = fullfile(subject.path,'anat');
if exist(pth,'dir')
    a = spm_select('List',pth,...
        sprintf('^%s.*_([a-zA-Z0-9]+){1}\\.nii(\\.gz)?$',subject.name));
    if isempty(a), a = {}; else a = cellstr(a); end
    for i=1:numel(a)

        %-Anatomy imaging data file
        %------------------------------------------------------------------
        subject.anat(i).filename = a{i}; % or full path?
        labels = regexp(a{i},[...
            '^sub-[a-zA-Z0-9]+' ...            % sub-<participant_label>
            '(?<ses>_ses-[a-zA-Z0-9]+)?' ...   % ses-<label>
            '(?<acq>_acq-[a-zA-Z0-9]+)?' ...   % acq-<label>
            '(?<rec>_rec-[a-zA-Z0-9]+)?' ...   % rec-<label>
            '(?<fa>_fa-[a-zA-Z0-9]+)?' ...     % fa-<index>
            '(?<echo>_echo-[a-zA-Z0-9]+)?' ... % echo-<index>
            '(?<inv>_inv-[a-zA-Z0-9]+)?' ...   % inv-<index>
            '(?<run>_run-[a-zA-Z0-9]+)?' ...   % run-<index>
            '_(?<type>[a-zA-Z0-9]+){1}' ...    % type
            '\.nii(\.gz)?$'],'names');         % NIfTI file extension
        subject.anat(i).type = labels.type;
        subject.anat(i).ses  = regexprep(labels.ses,'^_[a-zA-Z0-9]+-','');
        subject.anat(i).acq  = regexprep(labels.acq,'^_[a-zA-Z0-9]+-','');
        subject.anat(i).rec  = regexprep(labels.rec,'^_[a-zA-Z0-9]+-','');
        subject.anat(i).fa   = regexprep(labels.fa,'^_[a-zA-Z0-9]+-','');
        subject.anat(i).echo = regexprep(labels.echo,'^_[a-zA-Z0-9]+-','');
        subject.anat(i).inv  = regexprep(labels.inv,'^_[a-zA-Z0-9]+-','');
        subject.anat(i).run  = regexprep(labels.run,'^_[a-zA-Z0-9]+-','');

    end
end

%--------------------------------------------------------------------------
%-Task imaging data
%--------------------------------------------------------------------------
pth = fullfile(subject.path,'func');
if exist(pth,'dir')
    
    %-Task imaging data file
    %----------------------------------------------------------------------
    f = spm_select('List',pth,...
        sprintf('^%s.*_task-.*_bold\\.nii(\\.gz)?$',subject.name));
    if isempty(f), f = {}; else f = cellstr(f); end
    for i=1:numel(f)

        %-Task imaging data file
        subject.func(i).filename = f{i}; % or full path?
        labels = regexp(f{i},[...
            '^sub-[a-zA-Z0-9]+' ...            % sub-<participant_label>
            '(?<ses>_ses-[a-zA-Z0-9]+)?' ...   % ses-<label>
            '_task-(?<task>[a-zA-Z0-9]+)?' ... % task-<task_label>
            '(?<acq>_acq-[a-zA-Z0-9]+)?' ...   % acq-<label>
            '(?<rec>_rec-[a-zA-Z0-9]+)?' ...   % rec-<label>
            '(?<fa>_fa-[a-zA-Z0-9]+)?' ...     % fa-<index>
            '(?<echo>_echo-[a-zA-Z0-9]+)?' ... % echo-<index>
            '(?<inv>_inv-[a-zA-Z0-9]+)?' ...   % inv-<index>
            '(?<run>_run-[a-zA-Z0-9]+)?' ...   % run-<index>
            '_(?<type>[a-zA-Z0-9]+){1}' ...    % type (bold|sbref)
            '\.nii(\.gz)?$'],'names');    % NIfTI file extension
        subject.func(i).type = labels.type;
        subject.func(i).ses  = regexprep(labels.ses,'^_[a-zA-Z0-9]+-','');
        subject.func(i).task = labels.task;
        subject.func(i).acq  = regexprep(labels.acq,'^_[a-zA-Z0-9]+-','');
        subject.func(i).rec  = regexprep(labels.rec,'^_[a-zA-Z0-9]+-','');
        subject.func(i).fa   = regexprep(labels.fa,'^_[a-zA-Z0-9]+-','');
        subject.func(i).echo = regexprep(labels.echo,'^_[a-zA-Z0-9]+-','');
        subject.func(i).inv  = regexprep(labels.inv,'^_[a-zA-Z0-9]+-','');
        subject.func(i).run  = regexprep(labels.run,'^_[a-zA-Z0-9]+-','');

    end
    
    %-Task events file
    %----------------------------------------------------------------------
    f = spm_select('List',pth,...
        sprintf('^%s.*_task-.*_events\\.tsv$',subject.name));
    if isempty(f), f = {}; else f = cellstr(f); end
    for i=1:numel(f)
        %-Task events file
        subject.func(end+1).filename = f{i}; % or full path?
        labels = regexp(f{i},[...
            '^sub-[a-zA-Z0-9]+' ...            % sub-<participant_label>
            '(?<ses>_ses-[a-zA-Z0-9]+)?' ...   % ses-<label>
            '_task-(?<task>[a-zA-Z0-9]+)?' ... % task-<task_label>
            '(?<acq>_acq-[a-zA-Z0-9]+)?' ...   % acq-<label>
            '(?<rec>_rec-[a-zA-Z0-9]+)?' ...   % rec-<label>
            '(?<fa>_fa-[a-zA-Z0-9]+)?' ...     % fa-<index>
            '(?<echo>_echo-[a-zA-Z0-9]+)?' ... % echo-<index>
            '(?<inv>_inv-[a-zA-Z0-9]+)?' ...   % inv-<index>
            '(?<run>_run-[a-zA-Z0-9]+)?' ...   % run-<index>
            '_events\.tsv$'],'names');         % NIfTI file extension
        subject.func(end).type = 'events';
        subject.func(end).ses  = regexprep(labels.ses,'^_[a-zA-Z0-9]+-','');
        subject.func(end).task = labels.task;
        subject.func(end).acq  = regexprep(labels.acq,'^_[a-zA-Z0-9]+-','');
        subject.func(end).rec  = regexprep(labels.rec,'^_[a-zA-Z0-9]+-','');
        subject.func(end).fa   = regexprep(labels.fa,'^_[a-zA-Z0-9]+-','');
        subject.func(end).echo = regexprep(labels.echo,'^_[a-zA-Z0-9]+-','');
        subject.func(end).inv  = regexprep(labels.inv,'^_[a-zA-Z0-9]+-','');
        subject.func(end).run  = regexprep(labels.run,'^_[a-zA-Z0-9]+-','');
        subject.func(end).meta = spm_load(fullfile(pth,f{i})); % ?
    end
        
    %-Physiological and other continuous recordings file
    %----------------------------------------------------------------------
    f = spm_select('List',pth,...
        sprintf('^%s.*_task-.*_(physio|stim)\\.tsv\.gz$',subject.name));
    % see also [_recording-<label>]
    if isempty(f), f = {}; else f = cellstr(f); end
    for i=1:numel(f)
        subject.func(end+1).filename = f{i}; % or full path?
        labels = regexp(f{i},[...
            '^sub-[a-zA-Z0-9]+' ...            % sub-<participant_label>
            '(?<ses>_ses-[a-zA-Z0-9]+)?' ...   % ses-<label>
            '_task-(?<task>[a-zA-Z0-9]+)?' ... % task-<task_label>
            '(?<acq>_acq-[a-zA-Z0-9]+)?' ...   % acq-<label>
            '(?<rec>_rec-[a-zA-Z0-9]+)?' ...   % rec-<label>
            '(?<fa>_fa-[a-zA-Z0-9]+)?' ...     % fa-<index>
            '(?<echo>_echo-[a-zA-Z0-9]+)?' ... % echo-<index>
            '(?<inv>_inv-[a-zA-Z0-9]+)?' ...   % inv-<index>
            '(?<run>_run-[a-zA-Z0-9]+)?' ...   % run-<index>
            '_(?<type>physio|stim){1}\.tsv\.gz$'],'names'); % NIfTI file extension
        subject.func(end).type = labels.type;
        subject.func(end).ses  = regexprep(labels.ses,'^_[a-zA-Z0-9]+-','');
        subject.func(end).task = labels.task;
        subject.func(end).acq  = regexprep(labels.acq,'^_[a-zA-Z0-9]+-','');
        subject.func(end).rec  = regexprep(labels.rec,'^_[a-zA-Z0-9]+-','');
        subject.func(end).fa   = regexprep(labels.fa,'^_[a-zA-Z0-9]+-','');
        subject.func(end).echo = regexprep(labels.echo,'^_[a-zA-Z0-9]+-','');
        subject.func(end).inv  = regexprep(labels.inv,'^_[a-zA-Z0-9]+-','');
        subject.func(end).run  = regexprep(labels.run,'^_[a-zA-Z0-9]+-','');
        
    end
end

%--------------------------------------------------------------------------
%-Fieldmap data
%--------------------------------------------------------------------------
pth = fullfile(subject.path,'fmap');
if exist(pth,'dir')
    f = spm_select('List',pth,...
        sprintf('^%s.*\\.nii(\\.gz)?$',subject.name));
    if isempty(f), f = {}; else f = cellstr(f); end
    j = 1;

    %-Phase difference image and at least one magnitude image
    %----------------------------------------------------------------------
    labels = regexp(f,[...
        '^sub-[a-zA-Z0-9]+' ...              % sub-<participant_label>
        '(?<ses>_ses-[a-zA-Z0-9]+)?' ...     % ses-<label>
        '(?<acq>_acq-[a-zA-Z0-9]+)?' ...     % acq-<label>
        '(?<run>_run-[a-zA-Z0-9]+)?' ...     % run-<index>
        '_phasediff\.nii(\.gz)?$'],'names'); % NIfTI file extension
    if any(~cellfun(@isempty,labels))
        idx = find(~cellfun(@isempty,labels));
        for i=1:numel(idx)
            fb = spm_file(spm_file(f{idx(i)},'basename'),'basename');
            metafile = fullfile(pth,spm_file(fb,'ext','json'));
            subject.fmap(j).type = 'phasediff';
            subject.fmap(j).filename = f{idx(i)};
            subject.fmap(j).magnitude = {...
            	strrep(f{idx(i)},'_phasediff.nii','_magnitude1.nii'),...
            	strrep(f{idx(i)},'_phasediff.nii','_magnitude2.nii')}; % optional
            subject.fmap(j).ses = regexprep(labels{idx(i)}.ses,'^_[a-zA-Z0-9]+-','');
            subject.fmap(j).acq = regexprep(labels{idx(i)}.acq,'^_[a-zA-Z0-9]+-','');
            subject.fmap(j).run = regexprep(labels{idx(i)}.run,'^_[a-zA-Z0-9]+-','');
            subject.fmap(j).meta = spm_jsonread(metafile);
            j = j + 1;
        end
    end

    %-Two phase images and two magnitude images
    %----------------------------------------------------------------------
    labels = regexp(f,[...
        '^sub-[a-zA-Z0-9]+' ...           % sub-<participant_label>
        '(?<ses>_ses-[a-zA-Z0-9]+)?' ...  % ses-<label>
        '(?<acq>_acq-[a-zA-Z0-9]+)?' ...  % acq-<label>
        '(?<run>_run-[a-zA-Z0-9]+)?' ...  % run-<index>
        '_phase1\.nii(\.gz)?$'],'names'); % NIfTI file extension
    if any(~cellfun(@isempty,labels))
        idx = find(~cellfun(@isempty,labels));
        for i=1:numel(idx)
            fb = spm_file(spm_file(f{idx(i)},'basename'),'basename');
            metafile = fullfile(pth,spm_file(fb,'ext','json'));
            subject.fmap(j).type = 'phase12';
            subject.fmap(j).filename = {...
            	f{idx(i)},...
            	strrep(f{idx(i)},'_phase1.nii','_phase2.nii')};
            subject.fmap(j).magnitude = {...
            	strrep(f{idx(i)},'_phase1.nii','_magnitude1.nii'),...
            	strrep(f{idx(i)},'_phase1.nii','_magnitude2.nii')};
            subject.fmap(j).ses = regexprep(labels{idx(i)}.ses,'^_[a-zA-Z0-9]+-','');
            subject.fmap(j).acq = regexprep(labels{idx(i)}.acq,'^_[a-zA-Z0-9]+-','');
            subject.fmap(j).run = regexprep(labels{idx(i)}.run,'^_[a-zA-Z0-9]+-','');
            subject.fmap(j).meta = {...
            	spm_jsonread(metafile),...
            	spm_jsonread(strrep(metafile,'_phase1.json','_phase2.json'))};
            j = j + 1;
        end
    end

    %-A single, real fieldmap image
    %----------------------------------------------------------------------
    labels = regexp(f,[...
        '^sub-[a-zA-Z0-9]+' ...             % sub-<participant_label>
        '(?<ses>_ses-[a-zA-Z0-9]+)?' ...    % ses-<label>
        '(?<acq>_acq-[a-zA-Z0-9]+)?' ...    % acq-<label>
        '(?<run>_run-[a-zA-Z0-9]+)?' ...    % run-<index>
        '_fieldmap\.nii(\.gz)?$'],'names'); % NIfTI file extension
    if any(~cellfun(@isempty,labels))
        idx = find(~cellfun(@isempty,labels));
        for i=1:numel(idx)
            fb = spm_file(spm_file(f{idx(i)},'basename'),'basename');
            metafile = fullfile(pth,spm_file(fb,'ext','json'));
            subject.fmap(j).type = 'fieldmap';
            subject.fmap(j).filename = f{idx(i)};
            subject.fmap(j).magnitude = strrep(f{idx(i)},'_fieldmap.nii','_magnitude.nii');
            subject.fmap(j).ses = regexprep(labels{idx(i)}.ses,'^_[a-zA-Z0-9]+-','');
            subject.fmap(j).acq = regexprep(labels{idx(i)}.acq,'^_[a-zA-Z0-9]+-','');
            subject.fmap(j).run = regexprep(labels{idx(i)}.run,'^_[a-zA-Z0-9]+-','');
            subject.fmap(j).meta = spm_jsonread(metafile);
            j = j + 1;
        end
    end

    %-Multiple phase encoded directions (topup)
    %----------------------------------------------------------------------
    labels = regexp(f,[...
        '^sub-[a-zA-Z0-9]+' ...          % sub-<participant_label>
        '(?<ses>_ses-[a-zA-Z0-9]+)?' ... % ses-<label>
        '(?<acq>_acq-[a-zA-Z0-9]+)?' ... % acq-<label>
        '_dir-(?<dir>[a-zA-Z0-9]+)?' ... % dir-<index>
        '(?<run>_run-[a-zA-Z0-9]+)?' ... % run-<index>
        '_epi\.nii(\.gz)?$'],'names');   % NIfTI file extension
    if any(~cellfun(@isempty,labels))
        idx = find(~cellfun(@isempty,labels));
        for i=1:numel(idx)
            fb = spm_file(spm_file(f{idx(i)},'basename'),'basename');
            metafile = fullfile(pth,spm_file(fb,'ext','json'));
            subject.fmap(j).type = 'epi';
            subject.fmap(j).filename = f{idx(i)};
            subject.fmap(j).ses = regexprep(labels{idx(i)}.ses,'^_[a-zA-Z0-9]+-','');
            subject.fmap(j).acq = regexprep(labels{idx(i)}.acq,'^_[a-zA-Z0-9]+-','');
            subject.fmap(j).dir = labels{idx(i)}.dir;
            subject.fmap(j).run = regexprep(labels{idx(i)}.run,'^_[a-zA-Z0-9]+-','');
            subject.fmap(j).meta = spm_jsonread(metafile);
            j = j + 1;
        end
    end
end

%--------------------------------------------------------------------------
%-MEG data
%--------------------------------------------------------------------------
pth = fullfile(subject.path,'meg');
if exist(pth,'dir')
    
    %-MEG data file
    %----------------------------------------------------------------------
    m = spm_select('List',pth,...
        sprintf('^%s.*_task-.*_meg\\..*[^json]$',subject.name));
    if isempty(m), m = {}; else m = cellstr(m); end
    for i=1:numel(m)

        %-MEG data file
        subject.meg(i).filename = m{i}; % or full path?
        labels = regexp(m{i},[...
            '^sub-[a-zA-Z0-9]+' ...            % sub-<participant_label>
            '(?<ses>_ses-[a-zA-Z0-9]+)?' ...   % ses-<label>
            '_task-(?<task>[a-zA-Z0-9]+)?' ... % task-<task_label>
            '(?<run>_run-[a-zA-Z0-9]+)?' ...   % run-<index>
            '(?<proc>_proc-[a-zA-Z0-9]+)?' ... % proc-<label>
            '_meg\..*$'],'names');
        subject.meg(i).type = 'meg';
        subject.meg(i).task = labels.task;
        subject.meg(i).ses  = regexprep(labels.ses,'^_[a-zA-Z0-9]+-','');
        subject.meg(i).run  = regexprep(labels.run,'^_[a-zA-Z0-9]+-','');
        subject.meg(i).proc = regexprep(labels.proc,'^_[a-zA-Z0-9]+-','');

    end
    
    %-MEG events file
    %----------------------------------------------------------------------
    m = spm_select('List',pth,...
        sprintf('^%s.*_task-.*_events\\.tsv$',subject.name));
    if isempty(m), m = {}; else m = cellstr(m); end
    for i=1:numel(m)
        subject.meg(end+1).filename = m{i}; % or full path?
        labels = regexp(m{i},[...
            '^sub-[a-zA-Z0-9]+' ...            % sub-<participant_label>
            '(?<ses>_ses-[a-zA-Z0-9]+)?' ...   % ses-<label>
            '_task-(?<task>[a-zA-Z0-9]+)?' ... % task-<task_label>
            '(?<run>_run-[a-zA-Z0-9]+)?' ...   % run-<index>
            '(?<proc>_proc-[a-zA-Z0-9]+)?' ... % proc-<label>
            '_events\.tsv$'],'names');
        subject.meg(end).type = 'events';
        subject.meg(end).task = labels.task;
        subject.meg(end).ses  = regexprep(labels.ses,'^_[a-zA-Z0-9]+-','');
        subject.meg(end).run  = regexprep(labels.run,'^_[a-zA-Z0-9]+-','');
        subject.meg(end).proc = regexprep(labels.proc,'^_[a-zA-Z0-9]+-','');
        subject.meg(end).meta = spm_load(fullfile(pth,m{i})); % ?
    end
        
    %-Channel description table
    %----------------------------------------------------------------------
    m = spm_select('List',pth,...
        sprintf('^%s.*_task-.*_channels\\.tsv$',subject.name));
    if isempty(m), m = {}; else m = cellstr(m); end
    for i=1:numel(m)
        subject.meg(end+1).filename = m{i}; % or full path?
        labels = regexp(m{i},[...
            '^sub-[a-zA-Z0-9]+' ...            % sub-<participant_label>
            '(?<ses>_ses-[a-zA-Z0-9]+)?' ...   % ses-<label>
            '_task-(?<task>[a-zA-Z0-9]+)?' ... % task-<task_label>
            '(?<run>_run-[a-zA-Z0-9]+)?' ...   % run-<index>
            '(?<proc>_proc-[a-zA-Z0-9]+)?' ... % proc-<label>
            '_channels\.tsv$'],'names');
        subject.meg(end).type = 'channels';
        subject.meg(end).task = labels.task;
        subject.meg(end).ses  = regexprep(labels.ses,'^_[a-zA-Z0-9]+-','');
        subject.meg(end).run  = regexprep(labels.run,'^_[a-zA-Z0-9]+-','');
        subject.meg(end).proc = regexprep(labels.proc,'^_[a-zA-Z0-9]+-','');
        subject.meg(end).meta = spm_load(fullfile(pth,m{i})); % ?
    end

    %-Session-specific files
    %----------------------------------------------------------------------
    m = spm_select('List',pth,...
        sprintf('^%s(_ses-[a-zA-Z0-9]+)?.*_(photo\\.jpg|fid\\.json|fidinfo\\.txt|headshape\\..*)$',subject.name));
    if isempty(m), m = {}; else m = cellstr(m); end
    for i=1:numel(m)
        subject.meg(end+1).filename = m{i}; % or full path?
        labels = regexp(m{i},[...
            '^sub-[a-zA-Z0-9]+' ...            % sub-<participant_label>
            '(?<ses>_ses-[a-zA-Z0-9]+)?' ...   % ses-<label>
            '(?<task>_task-[a-zA-Z0-9]+)?' ... % task-<task_label>
            '(?<run>_run-[a-zA-Z0-9]+)?' ...   % run-<index>
            '(?<proc>_proc-[a-zA-Z0-9]+)?' ... % proc-<label>
            '_(?<type>[a-zA-Z0-9]+)?\..*$'],'names');
        subject.meg(end).type = labels.type;
        subject.meg(end).task = regexprep(labels.task,'^_[a-zA-Z0-9]+-','');
        subject.meg(end).ses  = regexprep(labels.ses,'^_[a-zA-Z0-9]+-','');
        subject.meg(end).run  = regexprep(labels.run,'^_[a-zA-Z0-9]+-','');
        subject.meg(end).proc = regexprep(labels.proc,'^_[a-zA-Z0-9]+-','');
        subject.meg(end).meta = struct();
    end
end

%--------------------------------------------------------------------------
%-Behavioral experiment data
%--------------------------------------------------------------------------
pth = fullfile(subject.path,'beh');
if exist(pth,'dir')

    % There might be a [_ses-<session_label>]

    %-Event timing
    %----------------------------------------------------------------------
    f = spm_select('List',pth, ...
        sprintf('^%s_task-.*_events.tsv$',subject.name));
    if ~isempty(f)
        f = cellstr(f);
        for i=1:numel(f)
            subject.beh(1).events(i).filename = fullfile(pth,f{i});
            task = regexp(f{i},'.*task-([a-zA-Z0-9]+)_events\.tsv$','tokens');
            subject.beh(1).events(i).task = task{1}{1};
        end
    end

    %-Metadata
    %----------------------------------------------------------------------
    f = spm_select('List',pth, ...
        sprintf('^%s_task-.*_beh.json$',subject.name));
    if ~isempty(f)
        f = cellstr(f);
        for i=1:numel(f)
            subject.beh(1).meta(i).filename = fullfile(pth,f{i});
            task = regexp(f{i},'.*task-([a-zA-Z0-9]+)_beh.json$','tokens');
            subject.beh(1).meta(i).task = task{1}{1};
        end
    end

    %-Physiological recordings
    %----------------------------------------------------------------------
    f = spm_select('List',pth, ...
        sprintf('^%s_task-.*_physio.tsv.gz$',subject.name));
    if ~isempty(f)
        f = cellstr(f);
        for i=1:numel(f)
            subject.beh(1).physio(i).filename = fullfile(pth,f{i});
            task = regexp(f{i},'.*task-([a-zA-Z0-9]+)_physio\.tsv\.gz$','tokens');
            subject.beh(1).physio(i).task = task{1}{1};
            metafile = fullfile(pth,spm_file(spm_file(f{i},'basename'),'ext','json'));
            if exist(metafile,'file')
                subject.beh(1).physio(i).meta = spm_jsonread(metafile);
            else
                subject.beh(1).physio(i).meta = [];
            end
        end
    end

    %-Other continuous recordings
    %----------------------------------------------------------------------
    f = spm_select('List',pth, ...
        sprintf('^%s_task-.*_stim.tsv.gz$',subject.name));
    if ~isempty(f)
        f = cellstr(f);
        for i=1:numel(f)
            subject.beh(1).stim(i).filename = fullfile(pth,f{i});
            task = regexp(f{i},'.*task-([a-zA-Z0-9]+)_stim\.tsv\.gz$','tokens');
            subject.beh(1).stim(i).task = task{1}{1};
            metafile = fullfile(pth,spm_file(spm_file(f{i},'basename'),'ext','json'));
            if exist(metafile,'file')
                subject.beh(1).stim(i).meta = spm_jsonread(metafile);
            else
                subject.beh(1).stim(i).meta = [];
            end
        end
    end
end

%--------------------------------------------------------------------------
%-Diffusion imaging data
%--------------------------------------------------------------------------
pth = fullfile(subject.path,'dwi');
if exist(pth,'dir')
    f = spm_select('FPList',pth,...
        sprintf('^%s.*_([a-zA-Z0-9]+){1}\\.nii(\\.gz)?$',subject.name));
    if isempty(f), f = {}; else f = cellstr(f); end
    for i=1:numel(f)

        %-Diffusion imaging file
        %------------------------------------------------------------------
        p = parse_filename(f{i}, {'sub','ses','acq','run', 'bval','bvec'});
        subject.dwi = [subject.dwi p];

        %-bval file
        %------------------------------------------------------------------
        bvalfile = get_metadata(f{i},'^.*%s\\.bval$');
        if isfield(bvalfile,'filename')
            subject.dwi(end).bval = spm_load(bvalfile.filename); % ?
        end

        %-bvec file
        %------------------------------------------------------------------
        bvecfile = get_metadata(f{i},'^.*%s\\.bvec$');
        if isfield(bvalfile,'filename')
            subject.dwi(end).bvec = spm_load(bvecfile.filename); % ?
        end
        
    end
end


%--------------------------------------------------------------------------
%-Positron Emission Tomography imaging data
%--------------------------------------------------------------------------
pth = fullfile(subject.path,'pet');
if exist(pth,'dir')
    f = spm_select('List',pth,...
        sprintf('^%s.*_task-.*_pet\\.nii(\\.gz)?$',subject.name));
    if isempty(f), f = {}; else f = cellstr(f); end
    for i=1:numel(f)
        
        %-PET imaging file
        %------------------------------------------------------------------
        p = parse_filename(f{i}, {'sub','ses','task','acq','rec','run'});
        subject.pet = [subject.pet p];
        
    end
end


%==========================================================================
%-Perform a BIDS query
%==========================================================================
function result = BIDS_query(BIDS,query,varargin)
opts = parse_query(varargin);
switch query
%   case 'subjects'
%       result = regexprep(unique({BIDS.subjects.name}),'^[a-zA-Z0-9]+-','');
    case 'sessions'
        result = unique({BIDS.subjects.session});
        result = regexprep(result,'^[a-zA-Z0-9]+-','');
    case 'modalities'
        hasmod = arrayfun(@(y) structfun(@(x) isstruct(x) & ~isempty(x),y),...
            BIDS.subjects,'UniformOutput',false);
        hasmod = any([hasmod{:}],2);
        mods   = fieldnames(BIDS.subjects)';
        result = mods(hasmod);
    case {'subjects', 'tasks', 'runs', 'types', 'data', 'metadata'}
        %-Initialise output variable
        result = {};
        %-Filter according to subjects
        if any(ismember(opts(:,1),'sub'))
            subs = opts{ismember(opts(:,1),'sub'),2};
            opts(ismember(opts(:,1),'sub'),:) = [];
        else
            subs = unique({BIDS.subjects.name});
            subs = regexprep(subs,'^[a-zA-Z0-9]+-','');
        end
        %-Filter according to modality
        if any(ismember(opts(:,1),'modality'))
            mods = opts{ismember(opts(:,1),'modality'),2};
            opts(ismember(opts(:,1),'modality'),:) = [];
        else
            mods = BIDS_query(BIDS,'modalities');
        end
        %-Get optional target option for metadata query
        if strcmp(query,'metadata') && any(ismember(opts(:,1),'target'))
            target = opts{ismember(opts(:,1),'target'),2};
            opts(ismember(opts(:,1),'target'),:) = [];
            if iscellstr(target)
                target = substruct('.',target{1});
            end
        else
            target = [];
        end
        %-Perform query
        for i=1:numel(BIDS.subjects)                    
            if ~ismember(BIDS.subjects(i).name(5:end),subs), continue; end
            for j=1:numel(mods)
                d = BIDS.subjects(i).(mods{j});
                for k=1:numel(d)
                    sts = true;
                    for l=1:size(opts,1)
                        if ~isfield(d(k),opts{l,1}) || ~ismember(d(k).(opts{l,1}),opts{l,2})
                            sts = false;
                        end
                    end
                    switch query
                        case 'subjects'
                            if sts
                                result{end+1} = BIDS.subjects(i).name;
                            end
                        case 'data'
                            if sts && isfield(d(k),'filename')
                                result{end+1} = fullfile(BIDS.subjects(i).path,mods{j},d(k).filename);
                            end
                        case 'metadata'
                            if sts && isfield(d(k),'filename')
                                f = fullfile(BIDS.subjects(i).path,mods{j},d(k).filename);
                                result{end+1} = get_metadata(f);
                                if ~isempty(target)
                                    try
                                        result{end} = subsref(result{end},target);
                                    catch
                                        warning('Non-existent field for metadata.');
                                        result{end} = [];
                                    end
                                end
                            end
%                             if sts && isfield(d(k),'meta')
%                                 result{end+1} = d(k).meta;
%                             end
                        case 'runs'
                            if sts && isfield(d(k),'run')
                                result{end+1} = d(k).run;
                            end
                        case 'tasks'
                            if sts && isfield(d(k),'task')
                                result{end+1} = d(k).task;
                            end
                        case 'types'
                            if sts && isfield(d(k),'type')
                                result{end+1} = d(k).type;
                            end
                    end
                end
            end
        end
        %-Postprocessing output variable
        switch query
            case 'subjects'
                result = unique(result);
                result = regexprep(result,'^[a-zA-Z0-9]+-','');
            case 'data'
                result = result';
            case 'metadata'
                if numel(result) == 1
                    result = result{1};
                end
            case {'tasks','runs','types'}
                result = unique(result);
                result(cellfun('isempty',result)) = [];
        end
    otherwise
        error('Unable to perform BIDS query.');
end


%==========================================================================
%-Parse BIDS query
%==========================================================================
function query = parse_query(query)
if numel(query) == 1 && isstruct(query{1})
    query = [fieldnames(query{1}), struct2cell(query{1})];
else
    if mod(numel(query),2)
        error('Invalid input syntax.');
    end
    query = reshape(query,2,[])';
end
for i=1:size(query,1)
    if ischar(query{i,2})
        query{i,2} = cellstr(query{i,2});
    end
    for j=1:numel(query{i,2})
        if iscellstr(query{i,2})
            query{i,2}{j} = regexprep(query{i,2}{j},sprintf('^%s-',query{i,1}),'');
        end
    end
end


%==========================================================================
%-Parse filename
%==========================================================================
function p = parse_filename(filename,fields)
filename = spm_file(filename,'filename');
parts = strsplit(filename,'_');
p.filename = filename;
[p.type, p.ext] = strtok(parts{end},'.');
for i=1:numel(parts)-1
    d = strsplit(parts{i},'-');
    p.(d{1}) = d{2};
end
if nargin == 2
    for i=1:numel(fields)
        if ~isfield(p,fields{i})
            p.(fields{i}) = '';
        end
    end
    try
        p = orderfields(p,['filename','ext','type',fields]);
    catch
        warning('File %s does not match template.',filename);
    end
end


%==========================================================================
%-Get metadata
%==========================================================================
function meta = get_metadata(filename, pattern)
if nargin == 1, pattern = '^.*_%s\\.json$'; end
pth = fileparts(filename);
p = parse_filename(filename);

meta = struct();

if isfield(p,'ses') && ~isempty(p.ses)
    N = 4; % there is a session level in the hierarchy
else
    N = 3;
end
    
for n=1:N
    metafile = spm_select('FPList',pth, sprintf(pattern,p.type));
    if isempty(metafile), metafile = {}; else metafile = cellstr(metafile); end
    for i=1:numel(metafile)
        p2 = parse_filename(metafile{i});
        fn = setdiff(fieldnames(p2),{'filename','ext','type'});
        ismeta = true;
        for j=1:numel(fn)
            if ~isfield(p,fn{j}) || ~strcmp(p.(fn{j}),p2.(fn{j}))
                ismeta = false;
                break;
            end
        end
        if ismeta
            if strcmp(p2.ext,'.json')
                meta = update_metadata(meta,spm_jsonread(metafile{i}));
            else
                meta.filename = metafile{i};
            end
        end
    end
    pth = fullfile(pth,'..');
end


%==========================================================================
%-Inheritance principle
%==========================================================================
function s1 = update_metadata(s1,s2)
fn = fieldnames(s2);
for i=1:numel(fn)
    if ~isfield(s1,fn{i})
        s1.(fn{i}) = s2.(fn{i});
    end
end
