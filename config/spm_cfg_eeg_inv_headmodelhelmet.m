function headmodelhelmet = spm_cfg_eeg_inv_headmodelhelmet
% Configuration file for specifying the head model for source reconstruction
% This is for registration using new helmet design.
%__________________________________________________________________________
% Copyright (C) 2012-2016 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_inv_headmodelhelmet.m 6926 2016-11-09 22:13:19Z guillaume $


D = cfg_files;
D.tag = 'D';
D.name = 'M/EEG datasets';
D.filter = 'mat';
D.num = [1 Inf];
D.help = {'Select the M/EEG mat files.'};

val = cfg_entry;
val.tag = 'val';
val.name = 'Inversion index';
val.strtype = 'n';
val.help = {'Index of the cell in D.inv where the results will be stored.'};
val.val = {1};

comment = cfg_entry;
comment.tag = 'comment';
comment.name = 'Comment';
comment.strtype = 's';
comment.help = {'User-specified information about this inversion'};
comment.val = {''};

template = cfg_const;
template.tag = 'template';
template.name = 'Template';
template.val  = {1};
template.help = {''};

mri = cfg_files;
mri.tag = 'mri';
mri.name = 'Individual structural image';
mri.filter = 'image';
mri.ufilter = '.*';
mri.num     = [1 1];
mri.help = {'Select the subject''s structural image'};

cortex = cfg_files;
cortex.tag = 'cortex';
cortex.name = 'Custom cortical mesh';
cortex.filter = 'mesh';
cortex.ufilter = '.*';
cortex.num     = [0 1];
cortex.help = {'Select the subject''s cortical mesh. Leave empty for default'};
cortex.val = {{''}};

iskull = cfg_files;
iskull.tag = 'iskull';
iskull.name = 'Custom inner skull mesh';
iskull.filter = 'mesh';
iskull.ufilter = '.*';
iskull.num     = [0 1];
iskull.help = {'Select the subject''s inner skull mesh. Leave empty for default'};
iskull.val = {{''}};

oskull = cfg_files;
oskull.tag = 'oskull';
oskull.name = 'Custom outer skull mesh';
oskull.filter = 'mesh';
oskull.ufilter = '.*';
oskull.num     = [0 1];
oskull.help = {'Select the subject''s outer skull mesh. Leave empty for default'};
oskull.val = {{''}};

scalp = cfg_files;
scalp.tag = 'scalp';
scalp.name = 'Custom scalp mesh';
scalp.filter = 'mesh';
scalp.ufilter = '.*';
scalp.num     = [0 1];
scalp.help = {'Select the subject''s scalp mesh. Leave empty for default'};
scalp.val = {{''}};

custom = cfg_branch;
custom.tag = 'custom';
custom.name = 'Custom meshes';
custom.help = {'Provide custom individual meshes as GIfTI files'};
custom.val  = {mri, cortex, iskull, oskull, scalp};

meshes = cfg_choice;
meshes.tag = 'meshes';
meshes.name = 'Mesh source';
meshes.values = {template, mri, custom};
meshes.val = {template};
meshes.help = {''};

meshres = cfg_menu;
meshres.tag = 'meshres';
meshres.name = 'Mesh resolution';
meshres.help = {'Specify the resolution of the cortical mesh'};
meshres.labels = {'coarse', 'normal', 'fine'};
meshres.values = {1, 2, 3};
meshres.val = {2};

meshing = cfg_branch;
meshing.tag = 'meshing';
meshing.name = 'Meshes';
meshing.help = {'Create head meshes for building the head model'};
meshing.val  = {meshes, meshres};

fidname = cfg_entry;
fidname.tag = 'fidname';
fidname.name = 'M/EEG fiducial label';
fidname.strtype = 's';
fidname.help = {'Label of a fiducial point (as specified in the M/EEG dataset)'};

type = cfg_entry;
type.tag = 'type';
type.name = 'Type MRI coordinates';
type.strtype = 'r';
type.num = [1 3];
type.help = {'Type the coordinates (in MNI or native space depending on the MRI supplied) corresponding to the fiducial in the structural image.'};

fiducials_filename = fullfile(spm('dir'), 'EEGtemplates', 'fiducials.sfp');
fid = fopen(fiducials_filename ,'rt');
if fid == -1, error('Cannot open "%s".',fiducials_filename); end
fidtable = textscan(fid ,'%s %f %f %f');
fclose(fid);

select = cfg_menu;
select.tag = 'select';
select.name = 'Select from a list';
select.help = {'Select the corresponding fiducial point from a pre-specified list.'};
select.labels = fidtable{1}';
select.values = fidtable{1}';

specification = cfg_choice;
specification.tag = 'specification';
specification.name = 'How to specify?';
specification.values = {select, type};
specification.help = {''};

fiducial = cfg_branch;
fiducial.tag = 'fiducial';
fiducial.name = 'Fiducial';
fiducial.help = {'Specify fiducial for coregistration'};
fiducial.val  = {fidname, specification};

fiducials = cfg_repeat;
fiducials.tag = 'fiducials';
fiducials.name = 'Fiducials';
fiducials.help = {'Specify fiducials for coregistration (at least 3 fiducials need to be specified)'};
fiducials.num  = [3 Inf];
fiducials.values  = {fiducial};
fiducials.val = {fiducial fiducial fiducial};

useheadshape = cfg_menu;
useheadshape.tag = 'useheadshape';
useheadshape.name = 'Use headshape points?';
useheadshape.help = {'Use headshape points (if available)'};
useheadshape.labels = {'yes', 'no'};
useheadshape.values = {1, 0};
useheadshape.val = {0};

coregspecify = cfg_branch;
coregspecify.tag = 'coregspecify';
coregspecify.name = 'Specify coregistration parameters';
coregspecify.val = {fiducials, useheadshape};
coregspecify.help = {''};

% coregdefault = cfg_const;
% coregdefault.tag = 'coregdefault';
% coregdefault.name = 'Use pre-calulated native-MRI to dewar transform';
% coregdefault.help = {'Use pre-calulated native-MRI to dewar transform'};
% coregdefault.val  = {1};

coregdefault = cfg_files;
coregdefault.tag = 'coregdefault';
coregdefault.name = 'Custom native to dewar transfrom for subject''s coregdefault';
coregdefault.filter = 'helmet';
coregdefault.ufilter = '.*';
coregdefault.num     = [0 1];
coregdefault.help = {'Select the subject''s helmet to MEG dewar transform'};
coregdefault.val = {{'NOT DEFINED'}};

% coregerror = cfg_entry;
% coregerror.tag = 'coregerror';
% coregerror.name = 'Coregistration ERROR to add in mm';
% coregerror.strtype = 'r';
% coregerror.val = {0};
% coregerror.help = {'random coreg error to add to fiducuals-LEAVE AT ZERO UNLESS YOU ARE SURE ABOUT THIS'};
% 
% 

coregistration = cfg_choice;
coregistration.tag = 'coregistration';
coregistration.name = 'Coregistration';
coregistration.values = {coregspecify, coregdefault};
coregistration.val = {coregspecify};
coregistration.help = {'Coregistration'};

eeg = cfg_menu;
eeg.tag = 'eeg';
eeg.name = 'EEG head model';
eeg.help = {'Select the head model type to use for EEG (if present)'};
eeg.labels = {'EEG BEM', '3-Shell Sphere (experimental)'};
eeg.values = {'EEG BEM', '3-Shell Sphere (experimental)'};
eeg.val = {'EEG BEM'};

meg = cfg_menu;
meg.tag = 'meg';
meg.name = 'MEG head model';
meg.help = {'Select the head model type to use for MEG (if present)'};
meg.labels = {'Single Sphere', 'MEG Local Spheres', 'Single Shell'};
meg.values = {'Single Sphere', 'MEG Local Spheres', 'Single Shell'};
meg.val = {'Single Sphere'};

forward = cfg_branch;
forward.tag = 'forward';
forward.name = 'Forward model';
forward.val = {eeg, meg};
forward.help = {'Forward model'};

headmodelhelmet = cfg_exbranch;
headmodelhelmet.tag = 'headmodelhelmet';
headmodelhelmet.name = 'MEG helmet head model specification';
headmodelhelmet.val = {D, val, comment, meshing, coregistration, forward};
headmodelhelmet.help = {'Specify MEG head model for forward computation using helmet'};
headmodelhelmet.prog = @specify_headmodel;
headmodelhelmet.vout = @vout_specify_headmodel;
headmodelhelmet.modality = {'MEG'};

function  out = specify_headmodel(job)

out.D = {};

%- Loop over input datasets
%--------------------------------------------------------------------------

for i = 1:numel(job.D)
    
    D = spm_eeg_load(job.D{i});
    
    if ~isfield(D,'inv')
        val   = 1;
    elseif numel(D.inv)<job.val
        val   = numel(D.inv) + 1;
    else
        val   = job.val;
    end
    
    if  val ~= job.val
        error(sprintf('Cannot use the user-specified inversion index %d for dataset ', job.val, i));
    end
    
    D.val = val;
    
    %-Meshes
    %--------------------------------------------------------------------------
    if ~isfield(D,'inv')
        D.inv = {struct('mesh', [])};
    end
    
    D.inv{val}.date    = strvcat(date,datestr(now,15));
    D.inv{val}.comment = {job.comment};
    
    if isfield(job.meshing.meshes, 'template')
        sMRI = 1;
    elseif isfield(job.meshing.meshes, 'custom')
        sMRI = job.meshing.meshes.custom.mri{1};
    else
        sMRI = job.meshing.meshes.mri{1};
    end
    
    
    %sMRI='D:\Luzia\MQ0342.3_bakup2\wmsMQ0342-0003-00001-000001-01.img'
    D = spm_eeg_inv_mesh_ui(D, val, sMRI, job.meshing.meshres); %% the mesh on the custom (or template) mri
    % writes out new cortex mesh in native mri space
    spm_eeg_inv_checkmeshes(D);

    if isfield(job.meshing.meshes, 'custom')
        if ~isempty(job.meshing.meshes.custom.scalp{1})
            D.inv{val}.mesh.tess_scalp = job.meshing.meshes.custom.scalp{1};
        end
        
        if ~isempty(job.meshing.meshes.custom.oskull{1})
            D.inv{val}.mesh.tess_oskull = job.meshing.meshes.custom.oskull{1};
        end
        
        if ~isempty(job.meshing.meshes.custom.iskull{1})
            D.inv{val}.mesh.tess_iskull = job.meshing.meshes.custom.iskull{1};
        end
        

        if ~isempty(job.meshing.meshes.custom.cortex{1})
            
            D.inv{val}.mesh.tess_ctx = job.meshing.meshes.custom.cortex{1};
            
            defs.comp{1}.inv.comp{1}.def = {D.inv{val}.mesh.def};
            defs.comp{1}.inv.space = {D.inv{val}.mesh.sMRI};
            defs.out{1}.surf.surface = {D.inv{val}.mesh.tess_ctx};
            defs.out{1}.surf.savedir.savesrc = 1;
            
            out = spm_deformations(defs);
            
            D.inv{val}.mesh.tess_mni     = export(gifti(out.surf{1}), 'spm');
 
        end
    end
    
    
    %-Coregistration
    %--------------------------------------------------------------------------
    
    if isfield(job.coregistration, 'coregdefault')
        % register using the Troebinger helmet system
        H1=load(job.coregistration.coregdefault{1}); %%
        % a number of coordinate systems to reconcile here:
        % the helmet coregistration was made in dewar space (a ctf dewar based coordinate system)
        % the transformation from dewar space to the native MRI is given by dewDEFAULT2NATIVE
        % when not using fiducial coils one also has a default ctf head based coordinate system
        % to get from the default head centered coordinate system to dewar space use H1.defaultHead2MEGdewar
        %  to get from the default head centered coordinate system to the current (coil defined) head centered
        % coordinate system use defaultHead2currentHead
        
        dewDEFAULT2NATIVE=H1.MEGdewar2MRI'; % H1.MEG2MRI transforms from MEG default dewar space to native MRI
        
        
        try, % need to fix this- spm object changed since calibration
            nocoilpos=H1.Dnocoils.sensors('MEG').coilpos;
            nocoilfids=H1.Dnocoils.fiducials;
        catch
            nocoilpos=H1.Dnocoils.sensors.meg.coilpos;
            nocoilfids=H1.Dnocoils.fiducials;
            nocoillabels=H1.Dnocoils.sensors.meg.label;
        end;
        
        if max(max(abs(D.fiducials.fid.pnt-nocoilfids.fid.pnt)))>1e-3,
            error('Both fiducial and headcast info: reset with ''changeHeadPos -nominal'' on Acq machine'); %% NEED TO WORK ON THIS
            % changeHeadpos -nominal resets to default coil positions
            
        end;
        
        [c,ia,ib]=intersect(D.sensors('MEG').label,nocoillabels,'rows')
        defaultHead2currentHead=spm_eeg_inv_rigidreg(D.sensors('MEG').coilpos(ia,:)',nocoilpos(ib,:)');
        
        meegfid = D.fiducials;
        
        
        mrifid = [];
        
        mrifid = D.inv{val}.mesh.fid; %% fiducials in the native MRI space (obtained from inverse transform from standard space)

        

        megpts=meegfid.fid.pnt; %% fiducials in head (dewar/sensor) space
        % first convert these to points in default head centred space
        megpts_defaulthead=pinv(defaultHead2currentHead)*[megpts';ones(1,size(megpts,2))];
        % now convert these default-head centred points to dewar space
        megpts_dewar=H1.defaultHead2MEGdewar*megpts_defaulthead;
        % now from dewar points to native
        megpts_native=dewDEFAULT2NATIVE*megpts_dewar;

        
        % put all transforms into one big one: current head centred coordinates to native space
        currenthead2NATIVE=dewDEFAULT2NATIVE*H1.defaultHead2MEGdewar*pinv(defaultHead2currentHead);
        nativepts=currenthead2NATIVE*[megpts';ones(1,size(megpts,2))]
        
        mrilbl = meegfid.fid.label;
        mrifid.fid.pnt=nativepts(1:3,:)';
        mrifid.fid.label=mrilbl;
        
        D.inv{val}.mesh.fid=mrifid; %% set mri fid to be transformed MEG fid
%         if isfield(job,'coregerror');  %%
%             if job.coregerror>0,
%                 disp('ADDING COREG ERROR');
%             end;
%             meegfid.fid.pnt
%             
%             randn('seed',sum(100*clock));
%             
%             meegfid.fid.pnt=meegfid.fid.pnt+randn(size(meegfid.fid.pnt)).*job.coregerror;
%             meegfid.fid.pnt
%             
%             
%         else
%             coregerror=0;
%         end;
%         
        
        
        D = spm_eeg_inv_datareg_ui(D, D.val, meegfid, mrifid,0);
        
    else
        meegfid = D.fiducials;
        selection = spm_match_str(meegfid.fid.label, {job.coregistration.coregspecify.fiducial.fidname});
        meegfid.fid.pnt = meegfid.fid.pnt(selection, :);
        meegfid.fid.label = meegfid.fid.label(selection);
        
        mrifid = [];
        mrifid.pnt = D.inv{val}.mesh.fid.pnt;
        mrifid.fid.pnt = [];
        mrifid.fid.label = {job.coregistration.coregspecify.fiducial.fidname}';
        
        for j = 1:numel(job.coregistration.coregspecify.fiducial)
            if isfield(job.coregistration.coregspecify.fiducial(j).specification, 'select')
                lbl = job.coregistration.coregspecify.fiducial(j).specification.select;
                ind = strmatch(lbl, D.inv{val}.mesh.fid.fid.label);
                mrifid.fid.pnt(j, :) = D.inv{val}.mesh.fid.fid.pnt(ind, :);
            else
                mrifid.fid.pnt(j, :) = job.coregistration.coregspecify.fiducial(j).specification.type;
            end
        end
        
        D = spm_eeg_inv_datareg_ui(D, D.val, meegfid, mrifid, job.coregistration.coregspecify.useheadshape);
    end
    
    %-Compute forward model
    %----------------------------------------------------------------------
    D.inv{val}.forward = struct([]);
    
    for j = 1:numel(D.inv{val}.datareg)
        switch D.inv{val}.datareg(j).modality
            case 'EEG'
                D.inv{D.val}.forward(j).voltype = job.forward.eeg;
            case 'MEG'
                D.inv{D.val}.forward(j).voltype = job.forward.meg;
        end
    end
    
    D = spm_eeg_inv_forward(D);
    
    for j = 1:numel(D.inv{val}.forward)
        spm_eeg_inv_checkforward(D, D.val, j);
    end
    
    save(D);
    
    out.D{i, 1} = fullfile(D.path, D.fname);
end

function dep = vout_specify_headmodel(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'M/EEG dataset(s) with a forward model';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});
