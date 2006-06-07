function D = spm_eeg_inv_BSTfwdsol(S)

%=======================================================================
% FORMAT D = spm_eeg_inv_BSTparameters(S)
% 
% This function defines the required options for running the BrainStorm
% forward solution (bst_headmodeler.m)
%
% S		    - input struct
% (optional) fields of S:
% D			- filename of EEG/MEG mat-file
%
% Output:
% D			- EEG/MEG data struct (also written to files)
% 
% See also help lines in bst_headmodeler.m
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout & Christophe Phillips
% $Id: spm_eeg_inv_BSTfwdsol.m 547 2006-06-07 12:23:17Z john $

spm_defaults


try
    D = S; 
    clear S
catch
    D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
	D = spm_eeg_ldata(D);
end

val = length(D.inv);

UserDB(1).STUDIES = D.path;
UserDB(1).SUBJECTS = D.path;
UserDB(1).FILELIST = '';
StudyFile = '';
% Store this information in Matlab's preferences
setpref('BrainStorm','UserDataBase',UserDB);
setpref('BrainStorm','iUserDataBase',1);
   

OPTIONS = bst_headmodeler;

% Forward approach
OPTIONS.HeadModelFile = '';

if ~isfield(D,'modality')
    ModFlag = spm_input('Data type?','+1','EEG|MEG',[1 2]);
else
    if isempty(D.modality)
        ModFlag = spm_input('Data type?','+1','EEG|MEG',[1 2]);
    else
        ModFlag = D.modality;
    end
    
end

switch ModFlag
    case 1
        D.modality = 'EEG';
    case 2
        D.modality = 'MEG';
end   
    
if D.modality == 'MEG'
    MethodFlag = spm_input('Forward model','+1','1S|OvS',[1 2]);
    switch MethodFlag
        case 1
            OPTIONS.Method = 'meg_sphere';
        case 2
            OPTIONS.Method = 'meg_os';
    end
elseif D.modality == 'EEG'
    MethodFlag = spm_input('Forward model','+1','1S|3S|3SBerg|OvS',[1 2 3 4]);
    switch MethodFlag
        case 1
            OPTIONS.Method = 'eeg_sphere';
        case 2
            OPTIONS.Method = 'eeg_3sphere';
        case 3
            OPTIONS.Method = 'eeg_3sphereBerg';
        case 4
            OPTIONS.Method = 'eeg_os';
    end
end    


% I/O functions
[pth,nam,ext] = fileparts(D.inv{val}.mesh.sMRI);
OPTIONS.ImageGridFile = [nam '_BSTGainMatrix.bin'];
OPTIONS.ImageGridBlockSize = D.inv{val}.mesh.Ctx_Nv + 1;
OPTIONS.FileNamePrefix = '';
OPTIONS.Verbose = '0';


% Head Geometry (create tesselation file)
TessName = fullfile(pth, [nam '_tesselation.mat']);


load(D.inv{val}.mesh.tess_ctx);
% convert positions into m
if max(max(vert)) > 1 % non m
    if max(max(vert)) < 20 % cm
        vert = vert/100;
    elseif max(max(vert)) < 200 % mm
        vert = vert/1000;
    end
end        
Comment{1}      = 'Cortex Mesh';
Faces{1}        = face;
Vertices{1}     = vert';
Curvature{1}    = [];
VertConn{1}     = cell(D.inv{val}.mesh.Ctx_Nv,1);
for i = 1:D.inv{val}.mesh.Ctx_Nv
    [ii,jj]	= find(face == i);
    voi	= setdiff(unique(face(ii,:)),i);
    voi = reshape(voi,1,length(voi));
    VertConn{1}{i} = voi;
end

if ~isempty(D.inv{val}.mesh.tess_iskull)
    load(D.inv{val}.mesh.tess_iskull);
    % convert positions into m
    if max(max(vert)) > 1 % non m
        if max(max(vert)) < 20 % cm
            vert = vert/100;
        elseif max(max(vert)) < 200 % mm
            vert = vert/1000;
        end
    end        
    Comment{2}      = 'Inner Skull Mesh';
    Faces{2}        = face;
    Vertices{2}     = vert';
    Curvature{2}    = [];
    VertConn{2}     = [];
    indx = 3;
else
    indx = 2;
end

if ~isempty(D.inv{val}.mesh.tess_scalp)
    load(D.inv{val}.mesh.tess_scalp);
    % convert positions into m
    if max(max(vert)) > 1 % non m
        if max(max(vert)) < 20 % cm
            vert = vert/100;
        elseif max(max(vert)) < 200 % mm
            vert = vert/1000;
        end
    end        
    Comment{indx}      = 'Scalp Mesh';
    Faces{indx}        = face;
    Vertices{indx}     = vert';
    Curvature{indx}    = [];
    VertConn{indx}     = [];
    
    OPTIONS.Scalp.FileName = TessName;
    OPTIONS.Scalp.iGrid    = indx;
    
    % compute the best fitting sphere
    Iz = find( vert(:,3) > (max(vert(:,3)) + 2*min(vert(:,3)))/3 );
    if isempty(Iz)
        disp('Error in defining the Scalp Mesh');
        return
    end
    [Center,Radius]    = spm_eeg_inv_BestFitSph(vert(Iz,:));
    OPTIONS.HeadCenter = Center;
    OPTIONS.Radii      = [.88 .93 1]*Radius;
end
clear vert face norm;


D.inv{val}.forward.bst_tess = TessName;
if spm_matlab_version_chk('7.1')>=0
    save(D.inv{val}.forward.bst_tess,'-V6','Comment','Curvature','Faces','VertConn','Vertices');
else
    save(D.inv{val}.forward.bst_tess,'Comment','Curvature','Faces','VertConn','Vertices');
end

% Sensor Information
OPTIONS.ChannelFile = [nam '_BSTChannelFile.mat'];
OPTIONS.ChannelType = D.modality;
if ~isempty(D.inv{val}.datareg.sens_coreg)
    ChannelLoc = load(D.inv{val}.datareg.sens_coreg);
else
    ChannelLoc = load(spm_select(1, '.mat', 'Sensor coordinates in MRI native space'));
end
FldName             = fieldnames(ChannelLoc); 
sens                = getfield(ChannelLoc,FldName{1})';
if length(sens) ~= size(sens,2) & length(sens) > 3
    sens = sens';
end
% convert positions into m
if max(max(sens)) > 1 % non m
     if max(max(sens)) < 20 % cm
         sens = sens/100;
     elseif max(max(sens)) < 200 % mm
         sens = sens/1000;
     end
end        
% sensor orientations
if (D.modality == 'MEG')
    if isempty(D.inv{val}.datareg.sens_orient)
        D.inv{val}.datareg.sens_orient = spm_select(1, '.mat', 'Sensor orientations');
    end
    ChannelOrient = load(D.inv{val}.datareg.sens_orient);
    if ~isempty(ChannelOrient)
        FldName       = fieldnames(ChannelOrient); 
        orientation   = getfield(ChannelOrient,FldName{1})';
    end
end
for i = 1:length(sens)
    Channel(i) = struct('Loc',[],'Orient',[],'Comment','','Weight',[],'Type','','Name','');
    Channel(i).Loc = sens(:,i);
    if exist('orientation') == 1
        Channel(i).Orient = orientation(:,i);
    else
        Channel(i).Orient = [];
    end
    Channel(i).Comment = num2str(i);
    Channel(i).Weight = 1;
    Channel(i).Type = D.modality;
    Channel(i).Name = [D.modality ' ' num2str(i)];
end
    
save(fullfile(D.path,OPTIONS.ChannelFile),'Channel');
    
if strcmp(D.modality,'EEG')
    OPTIONS.ChannelLoc  = [];
end
clear sens

D.inv{val}.forward.bst_channel = fullfile(pth,OPTIONS.ChannelFile);

% Set OPTIONS.EEGRef here !


% Source Model
OPTIONS.SourceModel         = -1; % current dipole model
OPTIONS.Cortex.FileName     = TessName;
OPTIONS.Cortex.iGrid        = 1;
load(D.inv{val}.mesh.tess_ctx);
% convert positions into m
if max(max(vert)) > 1 % non m
    if max(max(vert)) < 20 % cm
        vert = vert/100;
    elseif max(max(vert)) < 200 % mm
        vert = vert/1000;
    end
end        
OPTIONS.GridLoc             = vert';
OPTIONS.GridOrient          = normal';
clear vert face norm
OPTIONS.ApplyGridOrient     = 1;    % 1 dipole per location (0 : 3 dipoles per location)
OPTIONS.VolumeSourceGrid    = 0;
OPTIONS.SourceLoc           = OPTIONS.GridLoc;
OPTIONS.SourceOrient        = OPTIONS.GridOrient;


% Forward computation
[G,OPTIONS] = bst_headmodeler(OPTIONS);
% PCA of the gain matrix
[Gnorm,VectP,ValP] = spm_eeg_inv_PCAgain(G);


% Remove big and useless matrices from the bst_otptions structure
OPTIONS.Channel       = [];
OPTIONS.ChannelLoc    = [];
OPTIONS.ChannelOrient = [];
OPTIONS.GridOrient    = [];
OPTIONS.SourceLoc     = [];
OPTIONS.SourceOrient  = [];
OPTIONS.GridLoc       = [];


% Savings
D.inv{val}.forward.bst_options = OPTIONS;
D.inv{val}.forward.gainmat     = fullfile(pth,[nam '_SPMgainmatrix_' num2str(val) '.mat']);
D.inv{val}.forward.pcagain     = fullfile(pth,[nam '_SPMgainmatrix_pca_' num2str(val) '.mat']);
if spm_matlab_version_chk('7.1') >=0
    save(D.inv{val}.forward.gainmat,'-V6','G');
else
    save(D.inv{val}.forward.gainmat,'G');
end
if spm_matlab_version_chk('7.1') >=0
    save(D.inv{val}.forward.pcagain,'-V6','Gnorm','VectP','ValP');
else
    save(D.inv{val}.forward.pcagain,'Gnorm','VectP','ValP');
end
clear G Gnorm VectP ValP

if spm_matlab_version_chk('7.1') >= 0
	save(fullfile(D.path, D.fname), '-V6', 'D');
else
	save(fullfile(D.path, D.fname), 'D');
end
