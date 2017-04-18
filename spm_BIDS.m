function BIDS = spm_BIDS(root)
% Parse directory structure formated according to the BIDS standard
% FORMAT BIDS = spm_BIDS(root)
% root   - directory formated according to BIDS
% BIDS   - structure containing the BIDS file layout
%__________________________________________________________________________
%
% BIDS (Brain Imaging Data Structure): http://bids.neuroimaging.io/
%   The brain imaging data structure, a format for organizing and
%   describing outputs of neuroimaging experiments.
%   K. J. Gorgolewski et al, Scientific Data, 2016.
%__________________________________________________________________________
% Copyright (C) 2016-2017 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_BIDS.m 7061 2017-04-18 16:57:36Z guillaume $


if ~nargin, root = pwd; end
root = spm_select('CPath',root);

%==========================================================================
%-BIDS structure
%==========================================================================

BIDS = struct(...
    'dir',root, ...             % BIDS directory
    'descr',struct([]), ...     % content of dataset_description.json
    'sessions',{{}},...         % cellstr of sessions
    'subjects',struct([]), ...  % structure array of subjects
    'tasks',{{}},...            % content of task-<task_label>_bold.json (to be merged with sub-/func/sub-_bold.json
    'scans',struct([]),...      % content of sub-<participant_label>_scans.tsv (should go within subjects)
    'sess',struct([]),...       % content of sub-participants_label>_sessions.tsv (should go within subjects)
    'participants',struct([])); % content of particpants.tsv

%==========================================================================
%-Validation of BIDS root directory
%==========================================================================
if isempty(BIDS.dir)
    error('A BIDS directory has to be specified.');
elseif ~exist(BIDS.dir,'dir')
    error('BIDS directory does not exist.');
elseif ~exist(fullfile(BIDS.dir,'dataset_description.json'),'file')
    error('BIDS directory not valid: missing dataset_description.json.');
end

%==========================================================================
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

%==========================================================================
%-Optional directories
%==========================================================================
% [code/]
% [derivatives/]
% [stimuli/]
% [sourcedata/]

%==========================================================================
%-Scans key file
%==========================================================================

% sub-<participant_label>/[ses-<session_label>/]
%     sub-<participant_label>_scans.tsv

%==========================================================================
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

%==========================================================================
%-Sessions file
%==========================================================================

% sub-<participant_label>/[ses-<session_label>/]
%      sub-<participant_label>[_ses-<session_label>]_sessions.tsv

%==========================================================================
%-Tasks
%==========================================================================
t = spm_select('FPList',BIDS.dir,...
    '^task-.*_(beh|bold|channels|physio|stim|meg)\.(json|tsv)$');
if isempty(t), t = {}; else t = cellstr(t); end
for i=1:numel(t)
    task = spm_file(t{i},'basename');    
    labels = regexp(task,[...
        '^task-(?<task>[a-zA-Z0-9]+)?' ... % task-<task_label>
        '(?<acq>_acq-[a-zA-Z0-9]+)?' ... % acq-<label>
        '(?<rec>_rec-[a-zA-Z0-9]+)?' ... % rec-<label>
        '(?<run>_run-[a-zA-Z0-9]+)?' ... % run-<index>
        '_(?<type>(beh|bold|channels|physio|stim|meg))$'],'names'); % type
    BIDS.tasks{i}.type = labels.type;
    BIDS.tasks{i}.task = labels.task;
    BIDS.tasks{i}.acq  = strrep(labels.acq,'_','');
    BIDS.tasks{i}.rec  = strrep(labels.rec,'_','');
    BIDS.tasks{i}.run  = strrep(labels.run,'_','');
    BIDS.tasks{i}.meta = spm_load(t{i});
end

%==========================================================================
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


%==========================================================================
%-Subject
%==========================================================================
function subject = parse_subject(p, subjname, sesname)

    subject.name = subjname;       % subject name ('sub-*')
    subject.path = fullfile(p,subjname,sesname); % full path to subject directory
    subject.session = sesname;     % session name ('' or 'ses-*')
    subject.anat = struct([]);     % anatomy imaging data
    subject.func = struct([]);     % task imaging data
    subject.fmap = struct([]);     % fieldmap data
    subject.beh = struct([]);      % behavioral experiment data
    subject.dwi = struct([]);      % diffusion imaging data
    subject.meg = struct([]);      % MEG data
    
    %----------------------------------------------------------------------
    %-Anatomy imaging data
    %----------------------------------------------------------------------
    pth = fullfile(subject.path,'anat');
    if exist(pth,'dir')
        a = spm_select('List',pth,...
            sprintf('^%s.*\\.nii(\\.gz)$',subject.name));
        if isempty(a), a = {}; else a = cellstr(a); end
        for i=1:numel(a)
            
            %-Anatomy imaging data file
            %--------------------------------------------------------------
            subject.anat(i).filename = a{i}; % or full path?
            labels = regexp(a{i},[...
                '^sub-[a-zA-Z0-9]+' ... % sub-<participant_label>
                '(?<ses>_ses-[a-zA-Z0-9]+)?' ... % ses-<label>
                '(?<acq>_acq-[a-zA-Z0-9]+)?' ... % acq-<label>
                '(?<rec>_rec-[a-zA-Z0-9]+)?' ... % rec-<label>
                '(?<run>_run-[a-zA-Z0-9]+)?' ... % run-<index>
                '_(?<type>[a-zA-Z0-9]+)?' ... % modality
                '\.nii(\.gz)?$'],'names'); % NIfTI file extension
            subject.anat(i).type = labels.type;
            subject.anat(i).ses = strrep(labels.ses,'_','');
            subject.anat(i).acq = strrep(labels.acq,'_','');
            subject.anat(i).rec = strrep(labels.rec,'_','');
            subject.anat(i).run = strrep(labels.run,'_','');
            
            %-Metadata file
            %--------------------------------------------------------------
            % (spm_file called twice to handle .nii and .nii.gz)
            metafile = fullfile(pth,spm_file(spm_file(a{i},'basename'),'ext','json'));
            if exist(metafile,'file')
                subject.anat(i).meta = spm_jsonread(metafile);
            else
                subject.anat(i).meta = '';
            end
        end
    end
    
    %----------------------------------------------------------------------
    %-Task imaging data
    %----------------------------------------------------------------------
    pth = fullfile(subject.path,'func');
    if exist(pth,'dir')
        f = spm_select('List',pth,...
            sprintf('^%s_task-.*_bold\\.nii(\\.gz)$',subject.name));
        if isempty(f), f = {}; else f = cellstr(f); end
        for i=1:numel(f)
            
            %-Task imaging data file
            %--------------------------------------------------------------
            % (spm_file called twice to handle .nii and .nii.gz)
            fb = spm_file(spm_file(f{i},'basename'),'basename');
            subject.func(i).filename = f{i}; % or full path?
            labels = regexp(f{i},[...
                '^sub-[a-zA-Z0-9]+' ... % sub-<participant_label>
                '(?<ses>_ses-[a-zA-Z0-9]+)?' ... % ses-<label>
                '_task-(?<task>[a-zA-Z0-9]+)?' ... % task-<task_label>
                '(?<acq>_acq-[a-zA-Z0-9]+)?' ... % acq-<label>
                '(?<rec>_rec-[a-zA-Z0-9]+)?' ... % rec-<label>
                '(?<run>_run-[a-zA-Z0-9]+)?' ... % run-<index>
                '_bold\.nii(\.gz)?$'],'names'); % NIfTI file extension
            subject.func(i).task = labels.task;
            subject.func(i).acq = strrep(labels.acq,'_','');
            subject.func(i).rec = strrep(labels.rec,'_','');
            subject.func(i).run = strrep(labels.run,'_','');
            
            %-Acquisition file
            %--------------------------------------------------------------
            metafile = fullfile(pth,spm_file(fb,'ext','json'));
            if exist(metafile,'file')
                subject.func(i).meta = spm_jsonread(metafile);
            else
                subject.func(i).meta = [];
            end
            
            %-Task events file
            %--------------------------------------------------------------
            eventsfile = fullfile(pth,spm_file(spm_file(fb(1:end-5),...
                'suffix','_events'),'ext','tsv'));
            if exist(eventsfile,'file')
                subject.func(i).events = spm_load(eventsfile);
            else
                subject.func(i).events = [];
            end
            
            %-Physiological and other continuous recordings file
            %--------------------------------------------------------------
            % see also [_recording-<label>]
            physiofile = fullfile(pth,spm_file(spm_file(fb(1:end-5),...
                'suffix','_physio'),'ext','tsv.gz'));
            if exist(physiofile,'file')
                subject.func(i).physio = spm_load(physiofile);
            else
                subject.func(i).physio = [];
            end
            % and metafile _physio.json
            stimfile = fullfile(pth,spm_file(spm_file(fb(1:end-5),...
                'suffix','_stim'),'ext','tsv.gz'));
            if exist(stimfile,'file')
                subject.func(i).stim = spm_load(stimfile);
            else
                subject.func(i).stim = [];
            end
            % and metafile _stim.json
        end
    end
    
    %----------------------------------------------------------------------
    %-Fieldmap data
    %----------------------------------------------------------------------
    pth = fullfile(subject.path,'fmap');
    if exist(pth,'dir')
        f = spm_select('List',pth,...
            sprintf('^%s.*\\.nii(\\.gz)$',subject.name));
        if isempty(f), f = {}; else f = cellstr(f); end
        j = 1;
        
        %-Phase difference image and at least one magnitude image
        %------------------------------------------------------------------
        labels = regexp(f,[...
            '^sub-[a-zA-Z0-9]+' ... % sub-<participant_label>
            '(?<ses>_ses-[a-zA-Z0-9]+)?' ... % ses-<label>
            '(?<acq>_acq-[a-zA-Z0-9]+)?' ... % acq-<label>
            '(?<run>_run-[a-zA-Z0-9]+)?' ... % run-<index>
            '_phasediff\.nii(\.gz)?$'],'names'); % NIfTI file extension
        if any(~cellfun(@isempty,labels))
            idx = find(~cellfun(@isempty,labels));
            for i=1:numel(idx)
                fb = spm_file(spm_file(f{idx(i)},'basename'),'basename');
                metafile = fullfile(pth,spm_file(fb,'ext','json'));
                subject.fmap{j} = struct(...
                    'phasediff',f{idx(i)},...
                    'magnitude1',strrep(f{idx(i)},'_phasediff.nii','_magnitude1.nii'),...
                    'magnitude2',strrep(f{idx(i)},'_phasediff.nii','_magnitude2.nii'),... % optional
                    'ses',strrep(labels{idx(i)}.ses,'_',''),...
                    'acq',strrep(labels{idx(i)}.acq,'_',''),...
                    'run',strrep(labels{idx(i)}.run,'_',''),...
                    'meta',spm_jsonread(metafile));
                j = j + 1;
            end
        end
        
        %-Two phase images and two magnitude images
        %------------------------------------------------------------------
        labels = regexp(f,[...
            '^sub-[a-zA-Z0-9]+' ... % sub-<participant_label>
            '(?<ses>_ses-[a-zA-Z0-9]+)?' ... % ses-<label>
            '(?<acq>_acq-[a-zA-Z0-9]+)?' ... % acq-<label>
            '(?<run>_run-[a-zA-Z0-9]+)?' ... % run-<index>
            '_phase1\.nii(\.gz)?$'],'names'); % NIfTI file extension
        if any(~cellfun(@isempty,labels))
            idx = find(~cellfun(@isempty,labels));
            for i=1:numel(idx)
                fb = spm_file(spm_file(f{idx(i)},'basename'),'basename');
                metafile = fullfile(pth,spm_file(fb,'ext','json'));
                subject.fmap{j} = struct(...
                    'phase1',f{idx(i)},...
                    'phase2',strrep(f{idx(i)},'_phase1.nii','_phase2.nii'),...
                    'magnitude1',strrep(f{idx(i)},'_phase1.nii','_magnitude1.nii'),...
                    'magnitude2',strrep(f{idx(i)},'_phase1.nii','_magnitude2.nii'),...
                    'ses',strrep(labels{idx(i)}.ses,'_',''),...
                    'acq',strrep(labels{idx(i)}.acq,'_',''),...
                    'run',strrep(labels{idx(i)}.run,'_',''),...
                    'meta1',spm_jsonread(metafile),...
                    'meta2',spm_jsonread(strrep(metafile,'_phase1.json','_phase2.json')));
                j = j + 1;
            end
        end
        
        %-A single, real fieldmap image
        %------------------------------------------------------------------
        labels = regexp(f,[...
            '^sub-[a-zA-Z0-9]+' ... % sub-<participant_label>
            '(?<ses>_ses-[a-zA-Z0-9]+)?' ... % ses-<label>
            '(?<acq>_acq-[a-zA-Z0-9]+)?' ... % acq-<label>
            '(?<run>_run-[a-zA-Z0-9]+)?' ... % run-<index>
            '_fieldmap\.nii(\.gz)?$'],'names'); % NIfTI file extension
        if any(~cellfun(@isempty,labels))
            idx = find(~cellfun(@isempty,labels));
            for i=1:numel(idx)
                fb = spm_file(spm_file(f{idx(i)},'basename'),'basename');
                metafile = fullfile(pth,spm_file(fb,'ext','json'));
                subject.fmap{j} = struct(...
                    'fieldmap',f{idx(i)},...
                    'magnitude',strrep(f{idx(i)},'_fieldmap.nii','_magnitude.nii'),...
                    'ses',strrep(labels{idx(i)}.ses,'_',''),...
                    'acq',strrep(labels{idx(i)}.acq,'_',''),...
                    'run',strrep(labels{idx(i)}.run,'_',''),...
                    'meta',spm_jsonread(metafile));
                j = j + 1;
            end
        end
        
        %-Multiple phase encoded directions (topup)
        %------------------------------------------------------------------
        labels = regexp(f,[...
            '^sub-[a-zA-Z0-9]+' ... % sub-<participant_label>
            '(?<ses>_ses-[a-zA-Z0-9]+)?' ... % ses-<label>
            '(?<acq>_acq-[a-zA-Z0-9]+)?' ... % acq-<label>
            '_dir-(?<dir>[a-zA-Z0-9]+)?' ... % dir-<index>
            '(?<run>_run-[a-zA-Z0-9]+)?' ... % run-<index>
            '_epi\.nii(\.gz)?$'],'names'); % NIfTI file extension
        if any(~cellfun(@isempty,labels))
            idx = find(~cellfun(@isempty,labels));
            for i=1:numel(idx)
                fb = spm_file(spm_file(f{idx(i)},'basename'),'basename');
                metafile = fullfile(pth,spm_file(fb,'ext','json'));
                subject.fmap{j} = struct(...
                    'epi',f{idx(i)},...
                    'ses',strrep(labels{idx(i)}.ses,'_',''),...
                    'acq',strrep(labels{idx(i)}.acq,'_',''),...
                    'dir',labels{idx(i)}.dir,...
                    'run',strrep(labels{idx(i)}.run,'_',''),...
                    'meta',spm_jsonread(metafile));
                j = j + 1;
            end
        end
    end
    
    %----------------------------------------------------------------------
    %-MEG data
    %----------------------------------------------------------------------
    pth = fullfile(subject.path,'meg');
    if exist(pth,'dir')
        m = spm_select('List',pth,...
            sprintf('^%s_task-.*_meg\\..*[^json]$',subject.name));
        if isempty(m), m = {}; else m = cellstr(m); end
        for i=1:numel(m)
            
            %-MEG data file
            %--------------------------------------------------------------
            fb = spm_file(m{i},'basename');
            subject.meg(1).data(i).filename = m{i}; % or full path?
            labels = regexp(m{i},[...
                '^sub-[a-zA-Z0-9]+' ... % sub-<participant_label>
                '(?<ses>_ses-[a-zA-Z0-9]+)?' ... % ses-<label>
                '_task-(?<task>[a-zA-Z0-9]+)?' ... % task-<task_label>
                '(?<run>_run-[a-zA-Z0-9]+)?' ... % run-<index>
                '(?<proc>_proc-[a-zA-Z0-9]+)?' ... % proc-<label>
                '_meg\..*$'],'names');
            subject.meg(1).data(i).task = labels.task;
            subject.meg(1).data(i).ses  = strrep(labels.ses,'_','');
            subject.meg(1).data(i).run  = strrep(labels.run,'_','');
            subject.meg(1).data(i).proc = strrep(labels.proc,'_','');
            
            %-Metadata file
            %--------------------------------------------------------------
            metafile = fullfile(pth,spm_file(fb,'ext','json'));
            if exist(metafile,'file')
                subject.meg(1).data(i).meta = spm_jsonread(metafile);
            else
                subject.meg(1).data(i).meta = [];
            end
            
            %-MEG events file
            %--------------------------------------------------------------
            eventsfile = fullfile(pth,spm_file(spm_file(fb(1:end-4),...
                'suffix','_events'),'ext','tsv'));
            if exist(eventsfile,'file')
                subject.meg(1).data(i).events = spm_load(eventsfile);
            else
                subject.meg(1).data(i).events = [];
            end
            
            %-Channel description table
            %--------------------------------------------------------------
            channelfile = fullfile(pth,spm_file(spm_file(fb(1:end-4),...
                'suffix','_channels'),'ext','tsv'));
            if exist(channelfile,'file')
                subject.meg(1).data(i).channels = spm_load(channelfile);
            else
                subject.meg(1).data(i).channels = [];
            end
        end
        
        %-Session-specific files
        %------------------------------------------------------------------
        m = spm_select('List',pth,...
            sprintf('^%s(_ses-[a-zA-Z0-9]+)?_(photo\\.jpg|fid\\.json|fidinfo\\.txt|headshape\\..*)$',subject.name));
        if isempty(m), m = {}; else m = cellstr(m); end
        for i=1:numel(m)
            subject.meg(1).landmarks(i).filename = m{i}; % or full path?
        end
    end
    
    %----------------------------------------------------------------------
    %-Behavioral experiment data
    %----------------------------------------------------------------------
    pth = fullfile(subject.path,'beh');
    if exist(pth,'dir')
        
        % There might be a [_ses-<session_label>]
        
        %-Event timing
        %------------------------------------------------------------------
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
        %------------------------------------------------------------------
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
        %------------------------------------------------------------------
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
        %------------------------------------------------------------------
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
    
    %----------------------------------------------------------------------
    %-Diffusion imaging data
    %----------------------------------------------------------------------
    pth = fullfile(subject.path,'dwi');
    if exist(pth,'dir')
        d = spm_select('List',pth,...
            sprintf('^%s.*_dwi\\.nii(\\.gz)$',subject.name));
        if isempty(d), d = {}; else d = cellstr(d); end
        for i=1:numel(d)
            
            %-Diffusion imaging file
            %--------------------------------------------------------------
            subject.dwi(i).filename = d{i}; % or full path?
            labels = regexp(d{i},[...
                '^sub-[a-zA-Z0-9]+' ... % sub-<participant_label>
                '(?<ses>_ses-[a-zA-Z0-9]+)?' ... % ses-<label>
                '(?<acq>_acq-[a-zA-Z0-9]+)?' ... % acq-<label>
                '(?<run>_run-[a-zA-Z0-9]+)?' ... % run-<index>
                '_dwi\.nii(\.gz)?$'],'names'); % NIfTI file extension
            subject.dwi(i).ses = strrep(labels.ses,'_','');
            subject.dwi(i).acq = strrep(labels.acq,'_','');
            subject.dwi(i).run = strrep(labels.run,'_','');
            
            %-Metadata file
            %--------------------------------------------------------------
            % (spm_file called twice to handle .nii and .nii.gz)
            metafile = spm_file(spm_file(a{i},'basename'),'ext','json');
            if exist(fullfile(pth,metafile),'file')
                subject.dwi(i).meta = spm_jsonread(metafile);
            else
                subject.dwi(i).meta = '';
            end
            
            %-bval file
            %--------------------------------------------------------------
            bvalfile = spm_file(spm_file(a{i},'basename'),'ext','bval');
            if exist(fullfile(pth,bvalfile),'file')
                subject.dwi(i).bval = spm_load(bvalfile);
            else
                subject.dwi(i).bval = [];
            end
            
            %-bvec file
            %--------------------------------------------------------------
            bvecfile = spm_file(spm_file(a{i},'basename'),'ext','bvec');
            if exist(fullfile(pth,bvecfile),'file')
                subject.dwi(i).bvec = spm_load(bvecfile);
            else
                subject.dwi(i).bvec = [];
            end
        end
    end


%-Inheritance principle
function s1 = update(s1,s2)
fn = fieldnames(s2);
for i=1:numel(fn)
    s1.fn{i} = s2.fn{i};
end
