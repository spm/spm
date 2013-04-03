function coregshift = spm_cfg_eeg_inv_coregshift
% configuration file for specifying the head model for source
% reconstruction. THis is to add deterministic or random displacements to
% simulate coregistration error. GRB
%_______________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_cfg_eeg_inv_coregshift.m 5381 2013-04-03 10:50:41Z gareth $

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



sdshift = cfg_entry;
sdshift.tag = 'sdshift';
sdshift.name = 'Standard deviation in x,y z in mm';
sdshift.strtype = 'r';
sdshift.num = [1 3];
sdshift.val = {[0 0 0]};
sdshift.help = {'The standard deviation (mm) of the random displacement to be added to all fiducials'};


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
coregshift.val = {D, val, meanshift, sdshift,pperror};
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
    mmshift=job.meanshift;
    
    if max(abs(job.sdshift)>0)|| abs((job.pperror)>0),
        disp('changing random seed and adding coreg error');
        randn('seed',sum(100*clock));
        mmshift=mmshift+job.sdshift.*randn(1,3);
        
        
    end;
    
    
    %[L0,D] = spm_eeg_lgainmat(D);
    meegfid.fid.pnt=meegfid.fid.pnt+repmat(mmshift,3,1)+randn(size(meegfid.fid.pnt)).*job.pperror;
       
%        % meegfid.fid.pnt=meegfid.fid.pnt+randn(size(meegfid.fid.pnt)).*job.coregerror;
% %     87.4190   36.3453   10.4075
% %   -25.1042   73.3485   12.6594
% %    18.1259  -64.6969   10.5064
%     %% NB just change the effective head model position rather than the actual fiducial locations
%     
     D = spm_eeg_inv_datareg_ui(D, D.val, meegfid, mrifid,0);
% %       1.0368    0.0926   -0.0473  -22.2709
% %    -0.0287    0.8593    0.4148  -37.9634
% %     0.0880   -0.4511    0.9678   17.7373
% %          0         0         0    1.0000
     D.inv{D.val}.gainmat=''; %% these will now be incorrect
      D = spm_eeg_inv_forward(D);
%      
      for j = 1:numel(D.inv{val}.forward)
          spm_eeg_inv_checkforward(D, D.val, j);
      end
% %     
%     [L1,D] = spm_eeg_lgainmat(D);
%   
    
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

