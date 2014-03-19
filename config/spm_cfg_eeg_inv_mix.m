function inv_mix = spm_cfg_eeg_inv_mix
% Configuration file for merging (using a new inversion) a number of
% imaging source inversion reconstructions
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_cfg_eeg_inv_mix.m 5924 2014-03-19 14:59:12Z gareth $

D = cfg_files;
D.tag = 'D';
D.name = 'M/EEG datasets';
D.filter = 'mat';
D.num = [1 Inf];
D.help = {'Select the M/EEG mat files or .mat files containing inverses'};

val = cfg_entry;
val.tag = 'val';
val.name = 'Inversion index';
val.strtype = 'n';
val.help = {'Index of the cell in D.inv (same for all files) where the forward model can be found and the results will be stored.'};
val.val = {1};

prefix = cfg_entry;
prefix.tag = 'prefix';
prefix.name = 'Merged file prefix';
prefix.strtype = 's';
prefix.val={'merged'};
prefix.help = {'Prefix for the new filename that will contain the merged inversion results'};





inv_mix = cfg_exbranch;
inv_mix.tag = 'inv_mix';
inv_mix.name = 'Merge source estimates from multiple inversions';
inv_mix.val = {D, val,prefix};
inv_mix.help = {'To merge different source level variance estimates based on the same data'};
inv_mix.prog = @run_inv_mix;
inv_mix.vout = @vout_inv_mix;
inv_mix.modality = {'MEG'};

function  out = run_inv_mix(job)


inverse = [];
if numel(job.D)<1,
    error('Need to add a number of files to combine');
end;


%% first compile multiple inversions
disp('Loading inversions');
allID=[];
allJ=[];
allqC=[];
allF=[];

surfdir=[];
for j=1:numel(job.D), %% move through files- assume that changing directory means surface(i.e. lead field also changes)
    [a1,b1,c1]=fileparts(job.D{j});
    surfdir=strvcat(surfdir,a1);
end;


gainfiles=[];
for j=1:numel(job.D), %% move through files- assume that changing directory means surface(i.e. lead field also changes)
    
    try
        spmfilename=job.D{j};
        D = spm_eeg_load(spmfilename);
        inv=D.inv{job.val};
    catch
        dum=load(job.D{j});
        spmfilename=dum.spmfilename;
        [a0,b1,c1]=fileparts(spmfilename);
        [a1,b0,c0]=fileparts(job.D{j});
        D=spm_eeg_load([a1 filesep b1 c1]);
        inv=dum.inv;
    end;
    
    allF(j)=inv.inverse.F;
    allqC(j,:)=inv.inverse.qC;
    
    allmesh{j}.M   = inv.mesh.tess_mni;
    allcortex{j}=inv.mesh.tess_ctx;
    %% transform back to temporal subspace of original data
    %% (as Us will be different for different surfaces- but just rotations)
    check_data(j,:)=inv.inverse.U{1}'*inv.inverse.Y(:,1);
    
    disp('Loading SPM gain matrices from surface directories');
    
    
    [L,D]=spm_eeg_lgainmat(D);
    gainfile=[deblank(surfdir(j,:)) filesep D.inv{D.val}.gainmat];
    gainfiles = strvcat(gainfiles,gainfile);
    allL1(j,:)=L(1,:);
    
    
end;


%% check original data was the same
if numel(job.D)>1,
    if max(std(check_data))>max(std(check_data'))/1e6,
        error('data is not the same for these files');
    end;
end;


%%NB OCCAMS RAZOR AT 3
useind=find(allF>max(allF)-3); %%
allqC=allqC(useind,:);
allmesh=allmesh(useind);
allcortex=allcortex(useind);
gainfiles=gainfiles(useind,:);
allL1=allL1(useind,:);
%% only take unique lead fields out
[dum,fileind,surfind]=unique(allL1,'rows');
ugainfiles=gainfiles(fileind,:);


[a1,b1,c1]=fileparts(spmfilename);
outdir=surfdir(1,1:max(strfind(surfdir(1,:),filesep))); %% go up one directory
outfilename=[outdir filesep job.prefix b1 c1];
disp(sprintf('Copying and renaming original SPM file to %s',outfilename));

D2=copy(D,outfilename);

if ~isfield(D2.inv{D2.val},'inverse'),
    disp('No inversion parameters in file, taking from last inversion');
    D2.inv{D2.val}.inverse=inv.inverse;
end;
    
clear dum inverse D

vert=[];
face=[];
cmap1=[];
cortexstr='';
for j=1:length(fileind),
    offset=uint32(size(vert,1).*ones(size(allmesh{fileind(j)}.M.face)));
    col=j*10;
    vert=[vert ;allmesh{fileind(j)}.M.vert];
    face=[face ;allmesh{fileind(j)}.M.face+offset];
    cmap1=[cmap1 ;repmat(col,size(allmesh{fileind(j)}.M.face,1),1)];
 
    cortexstr=[cortexstr sprintf('%s',allcortex{fileind(j)})];
    if j<length(fileind),
        cortexstr=[cortexstr ';'];
    end;
end;
figure;
h=trisurf(face,vert(:,1),vert(:,2),vert(:,3),cmap1)
set(h,'Linestyle','none');
alpha(0.1);

% UPDATE MNI MESH- MAY NEED TO UPDATE OTHERS TOO
 D2.inv{D2.val}.mesh.tess_mni.vert=vert;
 D2.inv{D2.val}.mesh.tess_mni.face=face;
 D2.inv{D2.val}.mesh.tess_ctx=cortexstr;

%D2.save;
clear vert face cmap1 offset


D2=spm_eeg_invert_classic_mix(D2,D2.val,allqC,surfind,ugainfiles);


disp(sprintf('Improvement in model evidence over best single solution %3.2f',D2.inv{D2.val}.inverse.F-max(allF)));



if ~iscell(D2)
    D2 = {D2};
end

for i = 1:numel(D2)
    save(D2{i});
end

out.D = {outfilename};

function dep = vout_inv_mix(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'M/EEG dataset(s) after imaging source reconstruction';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});

