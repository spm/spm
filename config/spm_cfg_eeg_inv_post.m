function inv_post = spm_cfg_eeg_inv_post
% Configuration file for taking a number of previous inversion results (maybe based on different data), smoothing and creating an approximate posterior
%
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_cfg_eeg_inv_post.m 6424 2015-04-24 14:33:33Z gareth $

D = cfg_files;
D.tag = 'D';
D.name = 'M/EEG datasets';
D.filter = 'mat';
D.num = [1 Inf];
D.help = {'Select the M/EEG file'};

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


postname = cfg_entry;
postname.tag = 'postname';
postname.name = 'Prefix for priors';
postname.strtype = 's';
postname.num = [1 Inf];
postname.val = {'priorset1'};
postname.help = {'Prefix for prior directory'};



inv_post = cfg_exbranch;
inv_post.tag = 'inv_post';
inv_post.name = 'Create approx posterior';
inv_post.val = {D, val,smth,postname};
inv_post.help = {'Use to combine posterior variance estimates from multiple sessions- possibly for use as prior in future session'};
inv_post.prog = @run_inv_post;
inv_post.vout = @vout_inv_post;
inv_post.modality = {'MEG'};

function  out = run_inv_post(job)


inverse = [];
if numel(job.D)~=1,
    error('Need a single dataset');
end;

D=spm_eeg_load(job.D{1});
val=job.val;

[a1,b1,c1]=fileparts(D.fname);
priordir=[D.path filesep job.postname '_' b1];

[priorfiles] = spm_select('FPListRec',priordir,'.*\.mat$');

Nfiles=size(priorfiles,1);
fprintf('Found %d prior files\n',Nfiles);
if Nfiles==0,
    error('No prior file found in directory: %s', priordir);
end;

mesh=D.inv{val}.mesh.tess_mni;

qsum=sparse(Nfiles,length(mesh.vert));

for j=1:Nfiles
    fprintf('Loading priorfile %d of %d \n',j,Nfiles);
    load(deblank(priorfiles(j,:)),'Qp','Qe','UL','F');
    
    [LCpL,Q,sumLCpL,QE,Cy,M,Cp,Cq,Lq]=spm_eeg_assemble_priors(UL,Qp,{Qe});
    allF(j)=F;
    
    qsum(j,:)=sqrt(diag(Cp)); %% sum standard deviation over source space rather than power
   
     
    
end;
sumq=sum(qsum);
sumq=sumq./sum(sumq); %% sums to unity

mesh.faces=mesh.face;mesh.vertices=mesh.vert;
mesh=rmfield(mesh,'face');
mesh=rmfield(mesh,'vert');
ssumq=spm_mesh_smooth(mesh,sumq',job.smth);

ssumq=ssumq/sum(ssumq);

mesh=spm_mesh_inflate(mesh);

figure;

trisurf(mesh.faces,mesh.vertices(:,1),mesh.vertices(:,2),mesh.vertices(:,3),ssumq);
colorbar;
title('smooth posterior');

%% write gifti

save(gifti(mesh), [priordir filesep 'post.surf.gii']);


fname = fullfile(priordir, ['post' '.gii']);

G = gifti;
G.private.metadata(1).name = 'SurfaceID';
G.private.metadata(1).value = ['post' '.surf.gii'];

G.cdata = full(ssumq);
G.cdata = G.cdata(:);

save(G, fname, 'ExternalFileBinary');

out.D = job.D;


function dep = vout_inv_post(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Smooth posterior from a number of imaging source reconstructions';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});

