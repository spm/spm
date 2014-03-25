function inv_post = spm_cfg_eeg_inv_post
% Configuration file for taking a number of previous inversion results (maybe based on different data), smoothing and creating an approximate posterior
%
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_cfg_eeg_inv_post.m 5926 2014-03-25 11:55:25Z gareth $

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

smth = cfg_entry;
smth.tag = 'smth';
smth.name = 'Smoothing';
smth.strtype = 'r';
smth.help = {'Smoothing of power over cortical surface'};
smth.val = {30};


prefix = cfg_entry;
prefix.tag = 'prefix';
prefix.name = 'Gifti file name';
prefix.strtype = 's';
prefix.val={'post1'};
prefix.help = {'name for gifti mesh that will contain the approximate posterior'};




inv_post = cfg_exbranch;
inv_post.tag = 'inv_post';
inv_post.name = 'Create approx posterior';
inv_post.val = {D, val,smth,prefix};
inv_post.help = {'Use to combine posterior variance estimates from multiple sessions- possibly for use as prior in future session'};
inv_post.prog = @run_inv_post;
inv_post.vout = @vout_inv_post;
inv_post.modality = {'MEG'};

function  out = run_inv_post(job)


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
    
    %%allF(j)=inv.inverse.F;
    qC=sqrt(inv.inverse.qC); %% sum standard deviation over source space rather than power
    ProbqC=qC./max(qC);
    if j==1,
        sumPqC=ProbqC;
        mesh= inv.mesh.tess_mni;
        
    else
        
        if max(max(mesh.vert-inv.mesh.tess_mni.vert))>0.01,
            error('all inversions need to lie on the same mesh');
        end;
        sumPqC=sumPqC+ProbqC;
    end;
    
end;

mesh.faces=mesh.face;mesh.vertices=mesh.vert;

sumPqC=spm_mesh_smooth(mesh,sumPqC,job.smth);
sumPqC=sumPqC./sum(sumPqC);
figure;
mesh=rmfield(mesh,'face');
mesh=rmfield(mesh,'vert');
mesh=spm_mesh_inflate(mesh);

figure;
trisurf(mesh.faces,mesh.vertices(:,1),mesh.vertices(:,2),mesh.vertices(:,3),sumPqC);
colorbar;
title('smooth posterior');

%% write gifti
save(gifti(mesh), [job.prefix '.surf.gii']);


fname = fullfile(pwd, [job.prefix '.gii']);

G = gifti;
G.private.metadata(1).name = 'SurfaceID';
G.private.metadata(1).value = [job.prefix '.surf.gii'];

G.cdata = full(sumPqC);
G.cdata = G.cdata(:);

save(G, fname, 'ExternalFileBinary');



out.D = {fname};

function dep = vout_inv_post(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Smooth posterior from a number of imaging source reconstructions';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});

