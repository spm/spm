function coregshift = spm_cfg_eeg_inv_coregshift
% configuration file for specifying the head model for source
% reconstruction. THis is to add deterministic or random displacements to
% simulate coregistration error. GRB
%_______________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_cfg_eeg_inv_coregshift.m 5422 2013-04-16 15:35:49Z gareth $

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



meanshift = cfg_entry;
meanshift.tag = 'meanshift';
meanshift.name = 'Displacement in x,y z in mm';
meanshift.strtype = 'r';
meanshift.num = [1 3];
meanshift.val = {[0 0 0]};
meanshift.help = {'The mean displacement (in meg space) of the headmodel in mm'};

meanangshift = cfg_entry;
meanangshift.tag = 'meanangshift';
meanangshift.name = 'The mean rotation (in meg space) in degrees';
meanangshift.strtype = 'r';
meanangshift.num = [1 3];
meanangshift.val = {[0 0 0]};
meanangshift.help = {'The mean rotation (in meg space) of the headmodel in degrees'};

sdshift = cfg_entry;
sdshift.tag = 'sdshift';
sdshift.name = 'Standard deviation in x,y z in mm';
sdshift.strtype = 'r';
sdshift.num = [1 3];
sdshift.val = {[0 0 0]};
sdshift.help = {'The standard deviation (mm) of the random displacement to be added to all fiducials'};

sdangshift = cfg_entry;
sdangshift.tag = 'sdangshift';
sdangshift.name = 'Standard deviation of rotation in degrees';
sdangshift.strtype = 'r';
sdangshift.num = [1 3];
sdangshift.val = {[0 0 0]};
sdangshift.help = {'The standard deviation (degrees) of the random rotation of the fiducials'};

pperror = cfg_entry;
pperror.tag = 'pperror';
pperror.name = 'Per point per dimension standard deviation of error (mm)';
pperror.strtype = 'r';
pperror.num = [1 1];
pperror.val = {[0]};
pperror.help = {'The standard deviation of the error (mm) on each fiducial in each dimension'};



coregshift = cfg_exbranch;
coregshift.tag = 'coregshift';
coregshift.name = 'Add head model error';
coregshift.val = {D, val, meanshift, sdshift,meanangshift,sdangshift,pperror};
coregshift.help = {'To simulate the effects of coregistration error'};
coregshift.prog = @specify_coregshift;
coregshift.vout = @vout_specify_coregshift;
coregshift.modality = {'MEG'};

function  out = specify_coregshift(job)

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
    if ~isfield(D,'inv'),
        error('no head model set up');
    end
    
    meegfid = D.fiducials;
    
    mrifid = D.inv{val}.mesh.fid; %% fiducials in the native MRI space (obtained from inverse transform from standard space)
    
    
    megpts=meegfid.fid.pnt; %% fiducials in head (dewar/sensor) space
    
    
   
    startpos=meegfid.fid.pnt;
    
    rshift=abs(job.sdshift)+abs(job.sdangshift)+abs(job.pperror);
    if rshift>0;
        disp('changing random seed and adding coreg error');
        randn('seed',sum(100*clock));
    end;
    
    P(1:3)= job.meanshift+randn(1,3).*job.sdshift; %%  TRanslation
    P(4:6)=(job.meanangshift+randn(1,3).*job.sdangshift).*pi/180;   %% rotation in radians
    
    [A] = spm_matrix(P); %% deterministic coreg error
    origfid=[meegfid.fid.pnt ones(size(meegfid.fid.pnt,1),1)]
    newfid=(A*origfid')'; %% translated and rotated
    newfid=newfid(:,1:3);
   
    meegfid.fid.pnt=newfid+randn(size(meegfid.fid.pnt)).*job.pperror; %% now adding random error to each point
    disp('Displaced fid');
       meegfid.fid.pnt

%     %% NB just change the effective head model position rather than the actual fiducial locations
%     
     D = spm_eeg_inv_datareg_ui(D, D.val, meegfid, mrifid,0);
     D.inv{D.val}.gainmat=''; %% these will now be incorrect
     
      D = spm_eeg_inv_forward(D);
     
%      
      for j = 1:numel(D.inv{val}.forward)
          spm_eeg_inv_checkforward(D, D.val, j);
      end
    
    save(D);
    
    out.D{i, 1} = fullfile(D.path, D.fname);
end

function dep = vout_specify_coregshift(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'M/EEG dataset(s) with a forward model';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});

