function eeg_pca_distort = spm_cfg_eeg_pca_distort
% Configuration file for creating distorted versions of subject anatomy
% Based on original antomical and predetermined 100 eigen component template space.
%__________________________________________________________________________

% Gareth Barnes
% Copyright (C) 2024 Imaging Neuroscience


eeg_pca_distort          = cfg_exbranch;
eeg_pca_distort.tag      = 'eeg_pca_distort';
eeg_pca_distort.name     = 'Distortion Trajectory';
eeg_pca_distort.val      = @eeg_pca_distort_cfg;
eeg_pca_distort.help     = {'To distort surface information in an M/EEG dataset based on typical normal variation'};
eeg_pca_distort.prog     = @specify_eeg_pca_distort;
eeg_pca_distort.vout     = @vout_specify_eeg_pca_distort;
eeg_pca_distort.modality = {'MEG'};


%==========================================================================
function varargout = eeg_pca_distort_cfg

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

PCAtemplate = cfg_files;
PCAtemplate.tag = 'PCAtemplate';
PCAtemplate.name = 'PCA Template location';
PCAtemplate.filter = 'nii';
PCAtemplate.val={''};
PCAtemplate.dir={'D:\matlab\JoseYael\Templates\'};
PCAtemplate.num = [1 1];
PCAtemplate.help = {'Select the generic PCA template file (Template_0.nii)'};

TemplateRedo       = cfg_menu;
TemplateRedo.tag    = 'TemplateRedo';
TemplateRedo.name   = 'Force recompute if exists';
TemplateRedo.val    = {'No'};
TemplateRedo.help   = {'Will force recompute template for this subject'}';
TemplateRedo.labels = {'Yes', 'No'}';
TemplateRedo.values = {'Yes', 'No'}';

DistortIndices = cfg_entry;
DistortIndices.tag = 'DistortIndices';
DistortIndices.name = 'Range of indices';
DistortIndices.strtype = 'i';
DistortIndices.num = [1 2];
DistortIndices.val = {[8 100]};
DistortIndices.help = {'The range (from, to) of PCA indices to distort'};

Zrange = cfg_entry;
Zrange.tag = 'Zrange';
Zrange.name = 'Range of Z';
Zrange.strtype = 'r';
Zrange.num = [1 2];
Zrange.val = {[-3 3]};
Zrange.help = {'The normal range (in Z) over which distortion happens (eg [-3,+3])'};

Npoints = cfg_entry;
Npoints.tag = 'Npoints';
Npoints.name = 'Number points';
Npoints.strtype = 'i';
Npoints.num = [1 1];
Npoints.val = {[17]};
Npoints.help = {'The number of points on the distortion trajectory (linearly spaced between z limits)'};

RandSeed = cfg_entry;
RandSeed.tag = 'RandSeed';
RandSeed.name = 'Random seed';
RandSeed.strtype = 'i';
RandSeed.num = [1 1];
RandSeed.val = {[1]};
RandSeed.help = {'The random seed that defines trajectory'};


[cfg,varargout{1}] = deal({D, val, PCAtemplate,TemplateRedo, DistortIndices, Npoints,Zrange,RandSeed});


%==========================================================================
function  out = specify_eeg_pca_distort(job)

out.D = {};

[templatedir]=fileparts(job.PCAtemplate{1}); %% get directory with PCA template info
Nb=job.Npoints; %% number of points on distortion trajectory
D = spm_eeg_load(job.D{1});
val=job.val;
%-Meshes
%----------------------------------------------------------------------
if ~isfield(D,'inv')
    error('no head model set up');
end

D.inv{val}.gainmat=''; %% these will now be incorrect

%D = spm_eeg_inv_forward(D);
nativeMRIname=D.inv{1}.mesh.sMRI
[a1,b1,c1]=fileparts(nativeMRIname)
latentfilename=[a1 filesep 'PCA' filesep 'latent_code.mat'];
if strcmp(job.TemplateRedo,'Yes'),
    if exist([a1 filesep 'PCA']),
        fprintf('\nRemoving exiting PCA directory for fresh start\n')
        rmdir([a1 filesep 'PCA'],'s');
    end;
end;
%% the space the MRI ends up in is not typical mni canonical space, but the template
%% space that was used to create the dictionary. Calling this Dtemplate
DtemplateMRIname=[a1 filesep 'PCA' filesep b1 c1];
if ~exist(latentfilename),
    fprintf('Did not find latent file, transforming')
    [z,output_folder]=spm_pca_get_transforms(nativeMRIname,templatedir);
else
    fprintf('\n Found existing latent file for this subject not transforming\n')
end;
mesh=[];
mesh.template = false;
mesh.sMRI     = DtemplateMRIname;
%% this doesn't mess with the Dtemplate MRI, but creates files (stored in mesh)
%% to get it back to classic MNI space
mesh          = spm_eeg_inv_spatnorm(mesh);



PCrange=[job.DistortIndices(1):job.DistortIndices(2)]; %% COMPONENTS THAT WILL BE DEFORMED
rng('default');


rng(job.RandSeed); %% same random deviations for everyone
sgnPC=sign(randn(1,length(PCrange))); %% this sets distortion trajectory

[output_folder,b1,c1]=fileparts(mesh.sMRI);

Mcan=D.inv{1}.mesh.tess_mni;
Dtemplatecortexname='mesh_cortex_template.gii';

Dt{1}=Dtemplatecortexname;
Tmesh         = spm_swarp(Mcan, mesh.def); %% put subject's canonical mesh (in mni space) into PCA template space
Dtfilename      = fullfile(output_folder, Dtemplatecortexname);
fprintf('\n Writing out subject mesh in template space: \n %s \n',Dtfilename)
save(gifti(Tmesh), Dtfilename);

if sum(job.Zrange)~=0,
    error('not set up for off-centre Z ranges at the moment');
    %% maybe in future.
end;
sp=abs(job.Zrange(1));

PCin=PCrange.*sgnPC;
fprintf('\nCreating a trajectory of %d brains from z=-%3.2f to %3.2f with seed %d\n',Nb,sp,sp,job.RandSeed)
[z0,z,archivos] = spm_pca_sample_brains(Nb,PCin , sp, output_folder, templatedir,Dt,job.RandSeed);		% Generate gii distorted brains

figure;
zvals=linspace(-sp,sp,Nb);
fullz=zeros(max(PCrange),Nb)
fullz(PCrange,:)=z0;
imagesc(sort(z0(1,:)),1:max(PCrange),fullz);

ylabel('Component');
xlabel('Trajectory');
dum=colorbar;
ylabel(dum,'z','FontSize',18,'Rotation',0)
set(gca,'Fontsize',18)


%figure;
Mnative=gifti(D.inv{1}.mesh.tess_ctx);
vnative			= Mnative.vertices; %% native vertices
%meshplot4spm(Mnative)

title('original')




figure;
for i1 = 1:Nb,
    g_smp	= gifti(archivos{i1});
    v3		= g_smp.vertices;
    
    d1=sqrt(dot((v3-vnative)',(v3-vnative)'));

    [dvals,dind]=sort(d1);
    
    
    h=trisurf(g_smp.faces,v3(:,1),v3(:,2),v3(:,3),d1);

    
    set(gca,'visible','off')
    
    dum=colorbar;
    ylabel(dum,'mm','FontSize',18,'Rotation',90);

    title(sprintf('z=%3.1f',zvals(i1)));
    set(gca,'Fontsize',18);
    
    v3=double(v3);
    
    err_dist_native2brian(i1)	= sqrt(3)*mean(rms(vnative-v3,2));	% Error vs non-distorted native space mesh (100 PCs)
    
    [n1,x1]=hist(d1,50);
    alld1(i1,:)=d1;
end; % for i1

figure;

x1=0:0.5:max(alld1(:));
[n1]=hist(alld1(1,:),x1);
[n3]=hist(alld1(round(Nb/4),:),x1);
[n4]=hist(alld1(round(Nb/2),:),x1);
h=plot(x1,n1,x1,n3,x1,n4);
set(h,'Linewidth',4);
set(gca,'Fontsize',18);
m1=mean(alld1(1,:)); m3=mean(alld1(round(Nb/4),:));m4=mean(alld1(round(Nb/2),:));
x1ind=min(find(x1>m1)); x3ind=min(find(x1>m3)); x4ind=min(find(x1>m4));
hold on;

legend1=sprintf('|z|=%3.1f, mean=%3.1f mm',abs(zvals(1)),m1);
legend3=sprintf('|z|=%3.1f, mean=%3.1f mm',abs(zvals(round(Nb/4))),m3);
legend4=sprintf('|z|=%3.1f, mean=%3.1f mm',zvals(round(Nb/2)),m4);
legend(strvcat(legend1,legend3,legend4))
xlabel('Deviation (mm)')
ylabel('Number vertices')

figure;
h3=plot(zvals,err_dist_native2brian,'x',0,0,'ko');
set(h3,'Linewidth',4);
xlabel('Trajectory')
ylabel('Mean deviation (mm)');
set(gca,'Fontsize',18);


save(D);

out.D{1} = fullfile(D.path, D.fname);
out.brain=D.inv{1}.mesh.tess_ctx;
out.brians=archivos;




%==========================================================================
function dep = vout_specify_eeg_pca_distort(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'M/EEG dataset(s) with a forward model';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});
