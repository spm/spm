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
% $Id$

spm_defaults

try
    D = S; 
    clear S
catch
    D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
	D = spm_eeg_ldata(D);
end

val = length(D.inv);

OPTIONS                 = bst_headmodeler;

% Forward approach
OPTIONS.HeadModelFile   = '';

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
OPTIONS.ImageGridFile = [nam '_BSTGainMatrix.mat'];
OPTIONS.ImageGridBlockSize = D.inv{val}.mesh.Ctx_Nv + 1;
OPTIONS.FileNamePrefix = '';
OPTIONS.Verbose = '0';


% Head Geometry (create tesselation file)
TessName = [nam '_tesselation.mat'];

load(D.inv{val}.mesh.tess_ctx);
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
    OPTIONS.HeadCenter = Center';
    OPTIONS.Radii      = [.88 .93 1]*Radius;
end
clear vert face norm;

D.inv{val}.forward.bst_tess = fullfile(pth,TessName);
save(fullfile(pth,TessName),'Comment','Curvature','Faces','VertConn','Vertices');


% Set OPTIONS.EEGRef here !


% Sensor Information
OPTIONS.ChannelFile = [nam '_BSTChannelFile.mat'];
OPTIONS.ChannelType = D.modality;
ChannelLoc          = load(D.inv{val}.datareg.sens_coreg);
FldName             = fieldnames(ChannelLoc); 
OPTIONS.ChannelLoc  = getfield(ChannelLoc,FldName{1})';

D.inv{val}.forward.bst_channel = fullfile(pth,OPTIONS.ChannelFile);


% Source Model
OPTIONS.SourceModel         = -1; % current dipole model
OPTIONS.Cortex.FileName     = TessName;
OPTIONS.Cortex.iGrid        = 1;
load(D.inv{val}.mesh.tess_ctx);
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


% Savings
D.inv{val}.forward.bst_options = OPTIONS;
D.inv{val}.forward.gainmat     = fullfile(pth,[nam '_SPMgainmatrix_' num2str(val) '.mat']);
D.inv{val}.forward.pcagain     = fullfile(pth,[nam '_SPMgainmatrix_pca_' num2str(val) '.mat']);
save(D.inv{val}.forward.gainmat,'G');
save(D.inv{val}.forward.pcagain,'Gnorm','VectP','ValP');
clear G Gnorm VectP ValP

save(fullfile(pth,D.fname),'D');
