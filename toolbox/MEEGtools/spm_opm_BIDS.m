function D = spm_opm_BIDS(S)
% Create a BIDS compatible file structure from an OPM M/EEG object
% FORMAT D = spm_opm_BIDS(S)
%   S               - input structure
%  fields of S:
%   S.D             - SPM MEEG object                 - Default: no Default
%   S.studyFolder   - path to study folder. If it     - Default: 'study'
%                     doesn't exist it will be created
%   S.subLabel      - unique subject label            - Default: no Default
%   S.sesLabel      - unique session label            - Default: no Default
%   S.taskLabel     - unique task label               - Default: no Default
%   S.runLabel      - unique run label                - Default: no Default
%   S.powerLineF    - power line frequency            - Default: 50
%   S.subAge        - subject Age                     - Default: 'n/a'
%   S.subGender     - Gender                          - Default: 'n/a'
%   S.subGroup      - subject Group                   - Default: 'n/a'
%   S.subHandedness - subject handedness              - Default: 'n/a'
% Output:
%  D           - The copied datafile is returned
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Tim Tierney
% $Id: spm_opm_BIDS.m 7414 2018-09-07 11:00:29Z spm $

%-Set default values
%--------------------------------------------------------------------------
errorMsg = 'an MEEG object must be supplied.';
if ~isfield(S, 'D'),            error(errorMsg); end
if ~isfield(S, 'studyFolder'),  S.studyFolder = 'study'; end
errorMsg = 'subLabel must be supplied.';
if ~isfield(S, 'subLabel'),     error(errorMsg); end
errorMsg = 'sesLabel must be supplied.';
if ~isfield(S, 'sesLabel'),     error(errorMsg); end
errorMsg = 'taskLabel must be supplied.';
if ~isfield(S, 'taskLabel'),     error(errorMsg); end
errorMsg= 'runLabel must be supplied.';
if ~isfield(S, 'runLabel'),     error(errorMsg); end
if ~isfield(S, 'powerLineF'),   S.powerLineF  = 50; end
if ~isfield(S, 'subAge'),       S.subAge = 'n/a'; end
if ~isfield(S, 'subGender'),    S.subGender = 'n/a'; end
if ~isfield(S, 'subGroup'),     S.subGroup = 'n/a'; end
if ~isfield(S, 'subHandedness'),S.subHandedness = 'n/a'; end


%-Study and Derivatives
%--------------------------------------------------------------------------

%location
pD = path(S.D);
stFold = fullfile(pD,S.studyFolder);

% study and derivatives Folder
studyExists= exist(stFold,'dir');
derivatives = fullfile(stFold,'derivatives');
if(~studyExists)
    mkdir(stFold);
    mkdir(derivatives);
end

%-Subject, session and meg folder
%--------------------------------------------------------------------------

subId = ['sub-',S.subLabel];
subFold = fullfile(pD,S.studyFolder,subId);
sessId= ['ses-',S.sesLabel];
sesFold = fullfile(subFold,sessId);
megfold = fullfile(sesFold,'meg');

% create Subject, session and MEG folder
mkdir(subFold);
mkdir(sesFold);
mkdir(megfold);

%-Check for anatomical image
%--------------------------------------------------------------------------
try
    anatPath =S.D.inv{1}.mesh.sMRI;
    anatFold  = fullfile(sesFold,'anat');
    mkdir(anatFold);
    [a,b,c]= fileparts(anatPath);
    
    if strmatch(c,'.img')
        h = spm_vol(anatPath);
        im = spm_read_vols(h);
        h2 = h;
        anatPrefx = [subId,'_',sessId,'_T1w','.nii'];
        anatFile = fullfile(anatFold,anatPrefx);
        h2.fname =anatFile;
        spm_write_vol(h2,im);
    else
        anatPrefx = [subId,'_',sessId,'_T1w',c];
        anatFile = fullfile(anatFold,[b,c]);
        copyfile(anatPath,anatFile);
    end
    
catch
    
end

%-Copy in MEG data to new folder
%--------------------------------------------------------------------------
megPrefx = [subId,'_',sessId,'_task-',S.taskLabel,'_run-',S.runLabel,'_meg.ds'];
megDataFold = fullfile(megfold,megPrefx);
mkdir(megDataFold);

%-Now all those json and tsv files
%--------------------------------------------------------------------------

% now copy the dataset in there
megFile = fullfile(megDataFold,[megPrefx(1:end-3),'.mat']);
D =  copy(S.D,megFile);

% channel struct
chanStruct = [];
chanStruct.name = chanlabels(D);
chanStruct.type = chantype(D);
chanStruct.units = units(D);
chanFilePath = fullfile(megfold,[megPrefx(1:end-6),'channels.tsv']);
spm_save(chanFilePath,chanStruct);

% meg Struct
megStruct = [];
megStruct.TaskName = S.taskLabel;
megStruct.SamplingFrequency = fsample(D);
megStruct.PowerLineFrequency = S.powerLineF ;
megStruct.SoftwareFilters = 'n/a' ;
megStruct.DewarPosition = 'n/a' ;
megStruct.DigitizedLandmarks = false ;
megStruct.DigitizedHeadPoints =  false;
megFilePath = fullfile(megfold,[megPrefx(1:end-6),'meg.json']);
spm_save(megFilePath,megStruct);


coStruct = [];
coStruct.MEGCoordinateSystem = 'Other';
coStruct.MEGCoordinateUnits = 'mm';
coStruct.MEGCoordinateSystemDescription = 'Nifti world coordinates(already corregistered with native MRI)';
coFilePath = fullfile(megfold,[megPrefx(1:end-6),'coordsystem.json']);
spm_save(coFilePath,coStruct);

% participants file
pFilePath = fullfile(stFold,'participants.tsv');
pFileExists = exist(pFilePath)==2;

if(pFileExists)
    pStruct= spm_load(pFilePath);
     if(~iscell(pStruct.Age))
            pStruct.Age=repmat({'n/a'},length(pStruct.Age),1);
        end
    
        if(~iscell(pStruct.gender))
            pStruct.gender=repmat({'n/a'},length(pStruct.gender),1);
        end
   
        
        if(~iscell(pStruct.group))
            pStruct.group=repmat({'n/a'},length(pStruct.group),1);
        end
       
        
        if(~iscell(pStruct.handedness))
            pStruct.handedness=repmat({'n/a'},length(pStruct.handedness),1);
        end
    
    IDalreadyExists = any(strcmp(pStruct.participant_id,subId));
    if(~IDalreadyExists)
        pStruct.participant_id = {pStruct.participant_id{:},subId}';
        pStruct.Age = {pStruct.Age{:},num2str(S.subAge)}';
        pStruct.gender = {pStruct.gender{:},S.subGender}';
        pStruct.group = {pStruct.group{:},S.subGroup}';
        pStruct.handedness = {pStruct.handedness{:},S.subHandedness}';
    end
else
    pStruct = [];
    pStruct.participant_id = {subId};
    pStruct.Age = {S.subAge};
    pStruct.gender = {S.subGender};
    pStruct.group = {S.subGroup};
    pStruct.handedness = {S.subHandedness};
end

spm_save(pFilePath,pStruct);

% dataset description file
dStruct = [];
dStruct.Name = S.studyFolder;
dStruct.BIDSVersion = '1.0.0-rc2';
dFilePath = fullfile(stFold,'dataset_description.json');
spm_save(dFilePath,dStruct);
