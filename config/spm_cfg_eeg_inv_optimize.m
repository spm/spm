function optsetup = spm_cfg_eeg_inv_optimize
%function optsetup = spm_cfg_eeg_inv_optimize
% configuration file to set up optimization routines for M/EEG source
% inversion
%_______________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes


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


postname = cfg_entry;
postname.tag = 'postname';
postname.name = 'Prefix for posterior';
postname.strtype = 's';
postname.num = [1 Inf];
postname.val = {'priorset1'};
postname.help = {'Prefix for prior directory'};

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
optsetup.val = {D, val,opttype,postname};
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
priordir=[D.path filesep job.postname '_' b1];

[priorfiles] = spm_select('FPListRec',priordir,'.*\.mat$');

Nfiles=size(priorfiles,1);
fprintf('Found %d prior files\n',Nfiles);

%% add in functional (from other modalities or experiment) hypotheses

Qe0=0;%% bounding ratio of noise to signal power
disp('NB NO min sensor noise level');  %% NO MIN SENSOR NOISE LEVEL
allF=zeros(Nfiles,1);
Qpall={};Qeall={};
Fmax=-Inf;


for j=1:Nfiles,
    load(deblank(priorfiles(j,:)),'Qp','Qe','UL','F');
    fprintf('Optimizing priorfile %d of %d on its own\n',j,Nfiles);
    if ~isempty(F),
        Fstart=F;
    else 
        Fstart=-Inf;
    end;
    
    for k=1:length(job.opttype),
        [F,M,Cq,Cp,Qe,Qp] = spm_eeg_invert_EBoptimise(AY,UL,job.opttype(k),Qp,Qe,Qe0);
    end;
    
    if F<Fstart,
        disp('Could not improve, returning to prior');
        load(deblank(priorfiles(j,:)),'Qp','Qe','UL','F');
    end;
        
     
    %% now try combining this solution with best so far
    if j>1, 
        disp('Running in addition to best prior so far');
        a=load(deblank(priorfiles(maxind,:)),'Qp','Qe','UL','F');
        Qpall=Qp; % posterior from current optimization
        for k=1:length(a.Qp), % 
            Qpall{end+1}=a.Qp{k}; % augment with posterior from best so far
        end;
        maxpriors=512; %% set a max limit
        if length(Qpall)>maxpriors,
            Qpall=Qpall{1:maxpriors};
            warning('Limiting priors');
        end;
        [F2,M2,Cq2,Cp2,Qe2,Qp2] = spm_eeg_invert_EBoptimise(AY,UL,job.opttype,Qpall,Qe,Qe0);
        if F2>Fmax, %% if the combined posterior is better then use this and resave in same file
            fprintf('F improvment of %3.2f, Combining posteriors from prior files %d and %d\n',F2-Fmax,maxind,j);
            F=F2;M=M2;Cq=Cq2;Qe=Qe2;Qp=Qp2;Cp=Cp2;  
        end;
    end;
    
    allF(j)=F;
    if F>Fmax,
         Fmax=F;
         maxind=j;
         figure;colormap('gray');
         spm_mip([diag(Cp)],mesh.vert,6);
         title(sprintf('best iter %d, F=%3.2f',j,F));colorbar;
         inverse.F=F;
         inverse.M=M;
         inverse.Cq=Cq;
         inverse.Cp=Cp; %% posterior source level
         inverse.Qe=Qe; %% posterior sensor level
         inverse.Qp=Qp;

     end;
     
    save(deblank(priorfiles(j,:)),'Qp','Qe','UL','F'); 
    
    
end; % for j
    
   


inverse.allF=allF;



%----------------------------------------------------------------------
% evaluate conditional expectation (of the sum over trials)
%----------------------------------------------------------------------
SSR   = 0;
SST   = 0;
J     = {};

UY=inverse.UY;

for j = 1:numel(UY),
    
    % trial-type specific source reconstruction
    %------------------------------------------------------------------
    J{j} = inverse.M*UY{j};
    
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

spm_eeg_inv_view_priors(priorfiles,mesh);

rmind=find(allF<max(allF)-3);
fprintf('Removing %d poorest prior files\n',length(rmind));
for j=1:length(rmind),
    priorfiles(rmind(j),:)
    delete(deblank(priorfiles(rmind(j),:)));
end;



out.D = job.D;


function dep = vout_opt_priors(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'M/EEG dataset(s) after imaging source reconstruction';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});


