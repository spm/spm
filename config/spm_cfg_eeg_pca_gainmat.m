function eeg_pca_gainmat = spm_cfg_eeg_pca_gainmat
% Configuration file for creating distorted versions of subject anatomy
% Based on original antomical and predetermined 100 eigen component template space.
%__________________________________________________________________________

% Gareth Barnes
% Copyright (C) 2024 Imaging Neuroscience


eeg_pca_gainmat          = cfg_exbranch;
eeg_pca_gainmat.tag      = 'eeg_pca_gainmat';
eeg_pca_gainmat.name     = 'Gain matrices for surfaces';
eeg_pca_gainmat.val      = @eeg_pca_gainmat_cfg;
eeg_pca_gainmat.help     = {'To compute new lead field/ gain matrices for multiple (disto'};
eeg_pca_gainmat.prog     = @specify_eeg_pca_gainmat;
eeg_pca_gainmat.vout     = @vout_specify_eeg_pca_gainmat;
eeg_pca_gainmat.modality = {'MEG'};


%==========================================================================
function varargout = eeg_pca_gainmat_cfg

persistent cfg
if ~isempty(cfg), varargout = {cfg}; return; end

D = cfg_files;
D.tag = 'D';
D.name = 'M/EEG datasets';
D.filter = 'mat';
D.val={''};
D.dir={'C:\Users\gbarnes\Documents\jimmydata\output\data\gb070167\'};
D.num = [1 1];
D.help = {'Select the M/EEG mat file'};

val = cfg_entry;
val.tag = 'val';
val.name = 'Inversion index';
val.strtype = 'n';
val.help = {'Index of the cell in D.inv where the results will be stored.'};
val.val = {1};


% surface_files = cfg_files;
% surface_files.tag = 'surface_files';
% surface_files.name = 'Surface files ';
% surface_files.filter = 'gii';
% surface_files.val={''};
% surface_files.dir={''};
% %surface_files.num = [Inf Inf];
% surface_files.help = {'Select the surface files (*.gii) for which to compute new lead fields'};

LFRedo       = cfg_menu;
LFRedo.tag    = 'LFRedo';
LFRedo.name   = 'Force recompute of lead fields';
LFRedo.val    = {'No'};
LFRedo.help   = {'Will force recompute lead fields even if they exist'};
LFRedo.labels = {'Yes', 'No'}';
LFRedo.values = {'Yes', 'No'}';

LFheadmodel       = cfg_menu;
LFheadmodel.tag    = 'LFheadmodel';
LFheadmodel.name   = 'Use which head models';
LFheadmodel.val    = {'Single Shell'};
LFheadmodel.help   = {'Select the head model to use'};
LFheadmodel.labels = {'Single Shell', 'Single Sphere'}';
LFheadmodel.values = {'Single Shell', 'Single Sphere'}';



[cfg,varargout{1}] = deal({D, val,LFheadmodel,LFRedo});


%==========================================================================
function  out = specify_eeg_pca_gainmat(job)


out.D = {};



val=job.val;
D = spm_eeg_load(job.D{val});

outdir=[D.path filesep 'outPCA' filesep];


Dc=D.copy(outdir); %% work with copy of the original

nativeMRIname=Dc.inv{1}.mesh.sMRI
[a1,b1,c1]=fileparts(nativeMRIname)

distortdir=[a1 filesep 'PCA' filesep 'Cerebros' filesep];
[surfacefiles] = spm_select('FPList',distortdir,'.*\.gii$');
Nfiles=size(surfacefiles,1);
fprintf('\n Found %d surfaces in %s \n',Nfiles,distortdir);
fprintf('\n Plus using original surface %s \n',Dc.inv{1}.mesh.tess_ctx);


%-Meshes
%----------------------------------------------------------------------
headmodel=job.LFheadmodel;


headmodeldir=[outdir filesep headmodel];
mkdir(headmodeldir); %% make a directory to contain headmodels and then surface leadfields
Dc.inv{val}.forward(1).voltype=headmodel;
Dc=spm_eeg_inv_forward(Dc);



for cerebro = 0:Nfiles, %% move over brians and brain
    fprintf('\n Creating lead fields for brian %d of %d \n',cerebro, Nfiles);



    if cerebro>0, %% keep file as it was on the first run

        %% replace subject's cortex and delete any ref to lead fields
        brianmesh = deblank(surfacefiles(cerebro,:));

    else
        %brianmesh=brainmesh; %% set to original mesh
        brianmesh=Dc.inv{val}.mesh.tess_ctx;
    end;
    [a1,gainmatname,c1]=fileparts(brianmesh);

    Dc.inv{val}.mesh.tess_ctx=brianmesh;

    Mnative=gifti(Dc.inv{val}.mesh.tess_ctx); %% native space
    Vn=Mnative.vertices; Nv=length(Vn);
    %% got this transform (native mri-> meg space) from  spm_eeg_inv_forward
    M    = Dc.inv{val}.datareg.fromMNI*Dc.inv{val}.mesh.Affine;
    M    = diag([1e-3 1e-3 1e-3 1])*M;
    Vctf=(M*[Mnative.vertices ones(Nv,1)]')';
    Dc.inv{val}.forward.mesh.vert=double(Vctf(:,1:3)); % topology stays the same, so no need to update faces
    fprintf('\nComputing lead fields')

    copygainmatfile=[headmodeldir filesep 'GainMat_' gainmatname '.mat'];
    if ~exist(copygainmatfile) || strcmp(job.LFRedo,'Yes'),
        Dc.inv{val}.gainmat = '';
        Dc.save;

        [G,Dc]	= spm_eeg_lgainmat(Dc);			% Generate/load lead field

        spm_copy([Dc.path filesep Dc.inv{val}.gainmat],copygainmatfile);
    end;
    allgainmat{cerebro+1}=copygainmatfile;
end;


%% reset original mesh in the copied file
Dc.inv{val}.forward.mesh.vert=Dc.inv{val}.forward.mesh.vert;
Dc.inv{val}.mesh.tess_ctx=D.inv{val}.mesh.tess_ctx;
Dc.inv{val}.gainmat = ''; %% remove LF reference
Dc.save;

out.Dc{1} = fullfile(Dc.path, Dc.fname); %% output copy of original file
out.gainmat=allgainmat;





%==========================================================================
function dep = vout_specify_eeg_pca_gainmat(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'M/EEG dataset(s) with a forward model';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});
