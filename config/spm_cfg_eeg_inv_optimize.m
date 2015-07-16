function optsetup = spm_cfg_eeg_inv_optimize
%function optsetup = spm_cfg_eeg_inv_optimize
% configuration file to set up optimization routines for M/EEG source
% inversion
%_______________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_cfg_eeg_inv_optimize.m 6499 2015-07-16 13:37:41Z gareth $


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
val.help = {'Index of the cell in D.inv where the forward model can be found and the results will be stored.'};
val.val = {1};


priorname = cfg_entry;
priorname.tag = 'priorname';
priorname.name = 'Prefix for priors';
priorname.strtype = 's';
priorname.num = [1 Inf];
priorname.val = {'priorset1'};
priorname.help = {'Prefix for prior directory'};

REMLopt = cfg_entry;
REMLopt.tag = 'REMLopt';
REMLopt.name = 'REML parameters';
REMLopt.strtype = 'r';
REMLopt.num = [1 2];
REMLopt.val = {[-4,16]};
REMLopt.help = {'Select REML parameters'};

ARDopt = cfg_entry;
ARDopt.tag = 'ARDopt';
ARDopt.name = 'Threshold for ARD hyperparameter';
ARDopt.strtype = 'r';
ARDopt.num = [1 1];
ARDopt.val = {[128]};
ARDopt.help = {'ARD threshold for pruning hyperparameters (relative to max)'};

GSopt = cfg_entry;
GSopt.tag = 'GSopt';
GSopt.name = 'Set number of greed search iterations';
GSopt.strtype = 'i';
GSopt.num = [1 1];
GSopt.val = {[16]};
GSopt.help = {'Number of greedy search iterations'};


opttype = cfg_repeat;
opttype.tag = 'opttype';
opttype.name = 'Optimize';
opttype.help = {'Specify the optimization scheme'};
opttype.num  = [1 Inf];
opttype.values  = {REMLopt,ARDopt,GSopt};
opttype.val  = {REMLopt};



optsetup = cfg_exbranch;
optsetup.tag = 'optsetup';
optsetup.name = 'Inversion optimization';
optsetup.val = {D,val,priorname,opttype};
optsetup.help = {'Set optimization scheme for source reconstruction'};
optsetup.prog = @opt_priors;
optsetup.vout = @vout_opt_priors;

function  out = opt_priors(job)



D = spm_eeg_load(job.D{1});

%% get data specific terms sorted out
inverse=D.inv{job.val}.inverse;
AY=inverse.AY;

mesh=D.inv{job.val}.mesh.tess_mni;


[a1,b1,c1]=fileparts(D.fname);
priordir=[D.path filesep job.priorname '_' b1];


[priorfiles] = spm_select('FPListRec',priordir,'.*\.mat$');
%[postfiles] = spm_select('FPListRec',priordir,'^post.*\.mat$');

% USEPOST=1;
% if ~isempty(postfiles) && USEPOST,
%     fprintf('Using only posterior\n');
%     priorfiles=postfiles;
% end;
Nfiles=size(priorfiles,1);
fprintf('Found %d prior files\n',Nfiles);
if Nfiles==0,
    error('No prior file found in directory: %s', priordir);
end;

%% add in functional (from other modalities or experiment) hypotheses

Qe0=0;%% bounding ratio of noise to signal power
disp('NB NO min sensor noise level');  %% NO MIN SENSOR NOISE LEV

Qp_best=[];
Qe_best=[];
Fmax=-Inf;
F_aug=-Inf;
maxpriors=512;

for j=1:Nfiles, %% move through prior files
    %%% LOAD IN PRIOR FILE
    
    load(deblank(priorfiles(j,:)),'Qp','Qe','UL','F');
    fprintf('Optimizing priorfile %d of %d \n',j,Nfiles);
    %% ALSO CONSIDER AUGMENTED VERSION OF PRIORS IN FILE WITH BEST SO FAR
    Qp_aug=Qp; %% running with augmented Qp
    for k=1:length(Qp_best), %
        Qp_aug{length(Qp)+k}=Qp_best{k}; % augment with posterior from best so far
    end;
    Qe_aug=Qe; %% NB NOT AUGMENTING YET
    
    
    if length(Qp_aug)>maxpriors,
        Qp_aug=Qp_aug{1:maxpriors};
        warning('Limiting priors');
    end;
    %%%%%%%%%%%% NOW OPTIMIZE ORIGINAL AND AUGMENTED INDEPENDENTLY
    
    
    for k=1:length(job.opttype), %% move through optimization schemes
        %%% TERMS WHICH EVOLVE ARE Qe and Qp (M, Cq, Cp follow)
        priorcount(k)=length(Qp);
        
        [F,M,Cq,Cp,Qe,Qp] = spm_eeg_invert_EBoptimise(AY,UL,job.opttype(k),Qp,Qe,Qe0);
        
        if j>1, %% nothing to augment on 1st iteration
            
            [F_aug,M,Cq,Cp,Qe_aug,Qp_aug] = spm_eeg_invert_EBoptimise(AY,UL,job.opttype(k),Qp_aug,Qe_aug,Qe0);
        end;
    end;
    
    
    if F_aug>F,
        fprintf('Taking augmented set forward\n');
        Qp=Qp_aug;
        Qe=Qe_aug;
        F=F_aug;
    end;
    
    if F>Fmax, %% keep a record of these if they are best
        Qp_best=Qp;
        Qe_best=Qe;
        Fmax=F;
    end;
    
    Qp=Qp_best;
    Qe=Qe_best;
    save(deblank(priorfiles(j,:)),'Qp','Qe','UL','F');
    allF(j)=F;
end; % for j


%% Now get M, Cp, Cq etc based on Qp and Qe

[LCpL,Q,sumLCpL,QE,Cy,M,Cp,Cq,Lq]=spm_eeg_assemble_priors(UL,Qp,{Qe});


inverse.F=Fmax;
inverse.M=M;
inverse.qC=Cq;
inverse.Cp=Cp; %% posterior source level
inverse.Qe=Qe; %% posterior sensor level
inverse.Qp=Qp;
inverse.Is=1:size(Cp,1); % %% temporary fix





%----------------------------------------------------------------------
% evaluate conditional expectation (of the sum over trials)
%----------------------------------------------------------------------
SSR   = 0;
SST   = 0;
J     = {};

UY=inverse.UY;



sourcevar=zeros(1,size(inverse.M,1));
for j = 1:numel(UY),
    
    % trial-type specific source reconstruction
    %------------------------------------------------------------------
    J{j} = inverse.M*UY{j};
    Jtime=J{j}*inverse.T'; %% J is Nvert* Nsamples
    sourcevar=sourcevar+var(Jtime'); %% ./inverse.qC';
    % sum of squares
    %------------------------------------------------------------------
    SSR  = SSR + sum(var((UY{j} - UL*J{j}))); %% changed variance calculation
    SST  = SST + sum(var( UY{j}));
    
end

inverse.J=J;


% accuracy; signal to noise (over sources)
%======================================================================
inverse.R2   = 100*(SST - SSR)/SST;
fprintf('Percent variance explained %.2f\n',full(inverse.R2));
D.inv{job.val}.inverse=inverse;



spm_eeg_invert_display(D);

rmind=find(allF<max(allF)-3);
fprintf('Removing %d poorest prior files\n',length(rmind));
for j=1:length(rmind),
    fprintf('Deleting %s\n',priorfiles(rmind(j),:));
    delete(deblank(priorfiles(rmind(j),:)));
end;






idnum=round(spm_data_id(sourcevar)*1000); %% get unique id for file
postfilename=[priordir filesep sprintf('post%d.mat',idnum)];
postgiftiname=[priordir filesep sprintf('post%d.gii',idnum)];
fprintf('Making posterior %s based on diagonal\n',postfilename);
Qp=[];
Qp{1}=sparse(diag(sourcevar));
Mmni=[];
Mmni.faces=D.inv{job.val}.mesh.tess_mni.face;
Mmni.vertices=D.inv{job.val}.mesh.tess_mni.vert;
Mmni.cdata=sourcevar';
Mmni=gifti(Mmni);


save(deblank(postfilename),'Qp','Qe','UL','Fmax','postgiftiname','-v7.3');
save(Mmni,postgiftiname);


figure;
spm_mip(sourcevar,mesh.vert',6);
title('Posterior');
colorbar;




D.save;
out.postname=postfilename;
out.postgiftiname=postgiftiname;
out.D = job.D;


function dep = vout_opt_priors(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'M/EEG dataset(s) after imaging source reconstruction';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});


