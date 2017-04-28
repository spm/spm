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
% $Id: spm_BIDS.m 7072 2017-04-28 16:52:36Z guillaume $


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
    'dir',root, ...             % BIDS directory
    'descr',struct([]), ...     % content of dataset_description.json
    'sessions',{{}},...         % cellstr of sessions
    'subjects',struct([]), ...  % structure array of subjects
    'tasks',struct([]),...      % content of task-<task_label>_bold.json (to be merged with sub-/func/sub-_bold.json
    'scans',struct([]),...      % content of sub-<participant_label>_scans.tsv (should go within subjects)
    'sess',struct([]),...       % content of sub-participants_label>_sessions.tsv (should go within subjects)
    'participants',struct([])); % content of participants.tsv

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
    BIDS.descr = spm_jsonread(fullfile(BIDS.dir,'dataset_description.json'));
    if ~isfield(BIDS.descr,'BIDSVersion') || ~isfield(BIDS.descr,'Name')
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

%-Tasks
%==========================================================================
t = spm_select('FPList',BIDS.dir,...
    '^task-.*_(beh|bold|channels|physio|stim|meg)\.(json|tsv)$');
if isempty(t), t = {}; else t = cellstr(t); end
for i=1:numel(t)
    task = spm_file(t{i},'basename');    
    labels = regexp(task,[...
        '^task-(?<task>[a-zA-Z0-9]+)?' ... % task-<task_label>
        '(?<acq>_acq-[a-zA-Z0-9]+)?' ...   % acq-<label>
        '(?<rec>_rec-[a-zA-Z0-9]+)?' ...   % rec-<label>
        '(?<run>_run-[a-zA-Z0-9]+)?' ...   % run-<index>
        '_(?<type>(beh|bold|channels|physio|stim|meg))$'],'names'); % type
    BIDS.tasks(i).type = labels.type;
    BIDS.tasks(i).task = labels.task;
    BIDS.tasks(i).acq  = regexprep(labels.acq,'^_[a-zA-Z0-9]+-','');
    BIDS.tasks(i).rec  = regexprep(labels.rec,'^_[a-zA-Z0-9]+-','');
    BIDS.tasks(i).run  = regexprep(labels.run,'^_[a-zA-Z0-9]+-','');
    BIDS.tasks(i).meta = spm_load(t{i});
end

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
            '_(?<type>[a-zA-Z0-9]+){1}' ...    % modality
            '\.nii(\.gz)?$'],'names');         % NIfTI file extension
        subject.anat(i).type = labels.type;
        subject.anat(i).ses  = regexprep(labels.ses,'^_[a-zA-Z0-9]+-','');
        subject.anat(i).acq  = regexprep(labels.acq,'^_[a-zA-Z0-9]+-','');
        subject.anat(i).rec  = regexprep(labels.rec,'^_[a-zA-Z0-9]+-','');
        subject.anat(i).fa   = regexprep(labels.fa,'^_[a-zA-Z0-9]+-','');
        subject.anat(i).echo = regexprep(labels.echo,'^_[a-zA-Z0-9]+-','');
        subject.anat(i).inv  = regexprep(labels.inv,'^_[a-zA-Z0-9]+-','');
        subject.anat(i).run  = regexprep(labels.run,'^_[a-zA-Z0-9]+-','');

        %-Metadata file
        %------------------------------------------------------------------
        % (spm_file called twice to handle .nii and .nii.gz)
        metafile = fullfile(pth,spm_file(spm_file(a{i},'basename'),'ext','json'));
        if exist(metafile,'file')
            subject.anat(i).meta = spm_jsonread(metafile);
        else
            subject.anat(i).meta = '';
        end
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
        % (spm_file called twice to handle .nii and .nii.gz)
        fb = spm_file(spm_file(f{i},'basename'),'basename');
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
            '_bold\.nii(\.gz)?$'],'names');    % NIfTI file extension
        subject.func(i).type = 'bold';
        subject.func(i).ses  = regexprep(labels.ses,'^_[a-zA-Z0-9]+-','');
        subject.func(i).task = labels.task;
        subject.func(i).acq  = regexprep(labels.acq,'^_[a-zA-Z0-9]+-','');
        subject.func(i).rec  = regexprep(labels.rec,'^_[a-zA-Z0-9]+-','');
        subject.func(i).fa   = regexprep(labels.fa,'^_[a-zA-Z0-9]+-','');
        subject.func(i).echo = regexprep(labels.echo,'^_[a-zA-Z0-9]+-','');
        subject.func(i).inv  = regexprep(labels.inv,'^_[a-zA-Z0-9]+-','');
        subject.func(i).run  = regexprep(labels.run,'^_[a-zA-Z0-9]+-','');

        %-Acquisition file
        metafile = fullfile(pth,spm_file(fb,'ext','json'));
        if exist(metafile,'file')
            subject.func(i).meta = spm_jsonread(metafile);
        else
            subject.func(i).meta = [];
        end
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
        
    %-Physiological and other continuous recordings file (TODO)
    %----------------------------------------------------------------------
%     f = spm_select('List',pth,...
%         sprintf('^%s.*_task-.*_physio\\.tsv\.gz$',subject.name));
%     if isempty(f), f = {}; else f = cellstr(f); end
%     for i=1:numel(f)
%         % see also [_recording-<label>]
%         physiofile = fullfile(pth,spm_file(spm_file(fb(1:end-5),...
%             'suffix','_physio'),'ext','tsv.gz'));
%         if exist(physiofile,'file')
%             subject.func(i).physio = spm_load(physiofile);
%         else
%             subject.func(i).physio = [];
%         end
%         % and metafile _physio.json
%         stimfile = fullfile(pth,spm_file(spm_file(fb(1:end-5),...
%             'suffix','_stim'),'ext','tsv.gz'));
%         if exist(stimfile,'file')
%             subject.func(i).stim = spm_load(stimfile);
%         else
%             subject.func(i).stim = [];
%         end
%         % and metafile _stim.json
%     end
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
        fb = spm_file(m{i},'basename');
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

        %-Metadata file
        metafile = fullfile(pth,spm_file(fb,'ext','json'));
        if exist(metafile,'file')
            subject.meg(i).meta = spm_jsonread(metafile);
        else
            subject.meg(i).meta = [];
        end
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
        subject.meg(end).type = 'landmark';
        subject.meg(end).task = ''; % to fill in
        subject.meg(end).ses  = ''; % to fill in
        subject.meg(end).run  = ''; % to fill in
        subject.meg(end).proc = ''; % to fill in
        if strcmp(spm_file(m{i},'ext'),'json')
            subject.meg(end).meta = spm_jsonread(fullfile(pth,m{i}));
        else
            subject.meg(end).meta = [];
        end
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
    d = spm_select('List',pth,...
        sprintf('^%s.*_dwi\\.nii(\\.gz)?$',subject.name));
    if isempty(d), d = {}; else d = cellstr(d); end
    for i=1:numel(d)

        %-Diffusion imaging file
        %------------------------------------------------------------------
        subject.dwi(i).filename = d{i}; % or full path?
        labels = regexp(d{i},[...
            '^sub-[a-zA-Z0-9]+' ...          % sub-<participant_label>
            '(?<ses>_ses-[a-zA-Z0-9]+)?' ... % ses-<label>
            '(?<acq>_acq-[a-zA-Z0-9]+)?' ... % acq-<label>
            '(?<run>_run-[a-zA-Z0-9]+)?' ... % run-<index>
            '_dwi\.nii(\.gz)?$'],'names');   % NIfTI file extension
        subject.dwi(i).type = 'dwi';
        subject.dwi(i).ses = regexprep(labels.ses,'^_[a-zA-Z0-9]+-','');
        subject.dwi(i).acq = regexprep(labels.acq,'^_[a-zA-Z0-9]+-','');
        subject.dwi(i).run = regexprep(labels.run,'^_[a-zA-Z0-9]+-','');

        %-Metadata file
        %------------------------------------------------------------------
        % (spm_file called twice to handle .nii and .nii.gz)
        metafile = fullfile(pth,spm_file(spm_file(a{i},'basename'),'ext','json'));
        if exist(metafile,'file')
            subject.dwi(i).meta = spm_jsonread(metafile);
        else
            subject.dwi(i).meta = '';
        end

        %-bval file
        %------------------------------------------------------------------
        bvalfile = fullfile(pth,spm_file(spm_file(a{i},'basename'),'ext','bval'));
        if exist(fullfile(pth,bvalfile),'file')
            subject.dwi(i).bval = spm_load(bvalfile);
        else
            subject.dwi(i).bval = [];
        end

        %-bvec file
        %------------------------------------------------------------------
        bvecfile = fullfile(pth,spm_file(spm_file(a{i},'basename'),'ext','bvec'));
        if exist(fullfile(pth,bvecfile),'file')
            subject.dwi(i).bvec = spm_load(bvecfile);
        else
            subject.dwi(i).bvec = [];
        end
    end
end


%==========================================================================
%-Perform a BIDS query
%==========================================================================
function result = BIDS_query(BIDS,query,varargin)
opts = parse_query(varargin);
switch query
    case 'subjects'
        result = unique({BIDS.subjects.name});
        result = regexprep(result,'^[a-zA-Z0-9]+-','');
    case 'sessions'
        result = unique({BIDS.subjects.session});
        result = regexprep(result,'^[a-zA-Z0-9]+-','');
    case 'modalities'
        hasmod = arrayfun(@(y) structfun(@(x) isstruct(x) & ~isempty(x),y),...
            BIDS.subjects,'UniformOutput',false);
        hasmod = any([hasmod{:}],2);
        mods   = fieldnames(BIDS.subjects)';
        result = mods(hasmod);
    case 'types'
        result = {};
        mods = {'anat','func','fmap','dwi','meg'};
        for i=1:numel(BIDS.subjects)
            for j=1:numel(mods)
                for k=1:numel(BIDS.subjects(i).(mods{j}))
                    if isfield(BIDS.subjects(i).(mods{j})(k),'type')
                        result{end+1} = BIDS.subjects(i).(mods{j})(k).type;
                    end
                end
            end
        end
        result = unique(result);
    case 'tasks'
        result = unique({BIDS.tasks.task});
    case 'runs'
        result = {};
        if any(ismember(opts(:,1),'sub'))
            subs = opts{ismember(opts(:,1),'sub'),2};
            opts(ismember(opts(:,1),'sub'),:) = [];
        else
            subs = BIDS_query(BIDS,'subjects');
        end
        if any(ismember(opts(:,1),'modality'))
            mods = opts{ismember(opts(:,1),'modality'),2};
            opts(ismember(opts(:,1),'modality'),:) = [];
        else
            mods = BIDS_query(BIDS,'modalities');
        end
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
                    if sts && isfield(d(k),'run')
                    	result{end+1} = d(k).run;
                    end
                end
            end
        end
        result = unique(result);
        result(cellfun('isempty',result)) = [];
    case 'data'
        result = {};
        if any(ismember(opts(:,1),'sub'))
            subs = opts{ismember(opts(:,1),'sub'),2};
            opts(ismember(opts(:,1),'sub'),:) = [];
        else
            subs = BIDS_query(BIDS,'subjects');
        end
        if any(ismember(opts(:,1),'modality'))
            mods = opts{ismember(opts(:,1),'modality'),2};
            opts(ismember(opts(:,1),'modality'),:) = [];
        else
            mods = BIDS_query(BIDS,'modalities');
        end
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
                    if sts && isfield(d(k),'filename')
                    	result{end+1} = fullfile(BIDS.subjects(i).path,mods{j},d(k).filename);
                    end
                end
            end
        end
        result = result';
    case 'metadata'
        result = {};
        if any(ismember(opts(:,1),'sub'))
            subs = opts{ismember(opts(:,1),'sub'),2};
            opts(ismember(opts(:,1),'sub'),:) = [];
        else
            subs = BIDS_query(BIDS,'subjects');
        end
        if any(ismember(opts(:,1),'modality'))
            mods = opts{ismember(opts(:,1),'modality'),2};
            opts(ismember(opts(:,1),'modality'),:) = [];
        else
            mods = BIDS_query(BIDS,'modalities');
        end
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
                    if sts && isfield(d(k),'meta')
                    	result{end+1} = d(k).meta;
                    end
                end
            end
        end
    otherwise
        error('Unable to perform BIDS query.');
end


%-Parse BIDS query
function query = parse_query(query)
if numel(query) == 1 && isstruct(query{1})
    query = [fieldnames(query{1}), struct2cell(query{1})];
else
    query = reshape(query,2,[])';
end
for i=1:size(query,1)
    if ischar(query{i,2})
        query{i,2} = cellstr(query{i,2});
    end
    for j=1:numel(query{i,2})
        query{i,2}{j} = regexprep(query{i,2}{j},sprintf('^%s-',query{i,1}),'');
    end
end


%-Inheritance principle
function s1 = update(s1,s2)
fn = fieldnames(s2);
for i=1:numel(fn)
    s1.fn{i} = s2.fn{i};
end
