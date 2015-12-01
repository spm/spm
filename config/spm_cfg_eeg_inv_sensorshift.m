function sensorshift = spm_cfg_eeg_inv_sensorshift
% configuration file for tinkering with channel loations
%  THis is to add deterministic or random displacements to
% simulate coregistration error. GRB
%_______________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_cfg_eeg_inv_sensorshift.m 6618 2015-12-01 16:25:38Z spm $

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


movewhat = cfg_menu;
movewhat.tag = 'movewhat';
movewhat.name = 'Move all sensors together or sensors independently ?';
movewhat.help = {'Look at systematic errors in all sensor positioning or independent errors in individual sensors'};
movewhat.labels = {'all','independent'};
movewhat.values = {'all','independent'};
movewhat.val = {'independent'};


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

outprefix= cfg_entry;
outprefix.tag = 'outprefix';
outprefix.name = 'output prefix';
outprefix.strtype = 's';
outprefix.help = {'Prefix for new file with displaced sensors'};
outprefix.val = {'sens1'};%
% pperror = cfg_entry;
% pperror.tag = 'pperror';
% pperror.name = 'Per point per dimension standard deviation of error (mm)';
% pperror.strtype = 'r';
% pperror.num = [1 1];
% pperror.val = {[0]};
% pperror.help = {'The standard deviation of the error (mm) on each fiducial in each dimension'};



sensorshift = cfg_exbranch;
sensorshift.tag = 'sensorshift';
sensorshift.name = 'Add sensor error';
sensorshift.val = {D, val,movewhat, meanshift, sdshift,meanangshift,sdangshift,outprefix};
sensorshift.help = {'To simulate the effects of uncertainty on sensor position'};
sensorshift.prog = @specify_sensorshift;
sensorshift.vout = @vout_specify_sensorshift;
sensorshift.modality = {'MEG'};

function  out = specify_sensorshift(job)

out.D = {};

%- Loop over input datasets
%--------------------------------------------------------------------------

for i = 1:numel(job.D)
    
    D0 = spm_eeg_load(job.D{i});
    
    
    
    
    if ~isfield(D0,'inv')
        val   = 1;
    elseif numel(D0.inv)<job.val
        val   = numel(D0.inv) + 1;
    else
        val   = job.val;
    end
    
    if  val ~= job.val
        error(sprintf('Cannot use the user-specified inversion index %d for dataset ', job.val, i));
    end
    
    D0.val = val;
    
    
    %-Meshes
    %--------------------------------------------------------------------------
    if ~isfield(D0,'inv'),
        error('no head model set up');
    end;
    
    
    newfilename=[D0.path filesep job.outprefix D0.fname];
    D=D0.copy(newfilename);
    if isfield(D,'inv'),
        disp('Removing any previous inversions');
        D=rmfield(D,'inv');
    end;
    
    
    
    warning('OPERATING ON MEG SENSORS ONLY AS DEFAULT');
    
    
    sens1=D.sensors('MEG');
    
    chanind=D.indchantype('MEG');
    chanlabels=D.chanlabels(chanind);
    Nchans=length(chanind);
    origchanpos=D.sensors('MEG').chanpos;
    newchanpos=zeros(Nchans,3);
    newchanori=zeros(Nchans,3);
    
    switch job.movewhat
        case 'all',
            
            shift= job.meanshift+randn(1,3).*job.sdshift; %%  TRanslation
            rot=(job.meanangshift+randn(1,3).*job.sdangshift).*pi/180;   %% rotation in radians
            P=zeros(1,6);
            P(4:6)=rot;   %% rotation in radians
            [A] = spm_matrix(P); %% put rotation into matrix form
            for j=1:length(chanind),
                newchanpos(j,:)=sens1.chanpos(j,:)+shift; %%  TRanslation            
                newchanori(j,:)=sens1.chanori(j,:)*A(1:3,1:3)
            end; % for j
        case 'independent',
            for j=1:length(chanind),
                newchanpos(j,:)=sens1.chanpos(j,:)+job.meanshift+randn(1,3).*job.sdshift; %%  TRanslation
                P=zeros(1,6);
                P(4:6)=(job.meanangshift+randn(1,3).*job.sdangshift).*pi/180;   %% rotation in radians
                [A] = spm_matrix(P); %% put rotation into matrix form
                newchanori(j,:)=sens1.chanori(j,:)*A(1:3,1:3)
            end; % for j
    end;
    
    
    sens1.chanori=newchanori;
    sens1.coilori=newchanori;
    sens1.chanpos=newchanpos;
    sens1.coilpos=newchanpos;
    D=sensors(D,'MEG',sens1);
    
    D.save;
    
    
    out.D{i, 1} = fullfile(D.path, D.fname);
end

function dep = vout_specify_sensorshift(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'M/EEG dataset(s) with a forward model';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});

