function invert = spm_cfg_eeg_inv_invertiter
% Configuration file for configuring imaging source inversion reconstruction
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_cfg_eeg_inv_invertiter.m 6887 2016-09-20 08:01:36Z gareth $

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

outinv = cfg_entry;
outinv.tag = 'outinv';
outinv.name = 'Output prefix to save a copy of inv structure';
outinv.strtype = 's';
outinv.help = {'If this is supplied the original dataset will still be updated but new inverse structure with this name created'};utinv.val = {''};



all = cfg_const;
all.tag = 'all';
all.name = 'All';
all.val  = {1};

condlabel = cfg_entry;
condlabel.tag = 'condlabel';
condlabel.name = 'Condition label';
condlabel.strtype = 's';
condlabel.val = {''};

conditions = cfg_repeat;
conditions.tag = 'conditions';
conditions.name = 'Conditions';
conditions.help = {'Specify the labels of the conditions to be included in the inversion'};
conditions.num  = [1 Inf];
conditions.values  = {condlabel};
conditions.val = {condlabel};

whatconditions = cfg_choice;
whatconditions.tag = 'whatconditions';
whatconditions.name = 'What conditions to include?';
whatconditions.values = {all, conditions};
whatconditions.val = {all};

standard = cfg_const;
standard.tag = 'standard';
standard.name = 'Standard';
standard.help = {'Use default settings for the inversion'};
standard.val  = {1};

invfunc = cfg_menu;
invfunc.tag = 'invfunc';
invfunc.name = 'Inversion function';
invfunc.help = {'Current code allows multiple subjects and modalities; classic is for single subject single modality but has no additional scaling factors'};
invfunc.labels = {'Current','Classic'};
invfunc.values = {'Current','Classic'};
invfunc.val = {'Classic'};


invtype = cfg_menu;
invtype.tag = 'invtype';
invtype.name = 'Inversion type';
invtype.help = {'Select the desired inversion type'};
invtype.labels = {'GS', 'ARD', 'MSP (GS+ARD)' 'COH', 'IID','EBB','EBBsparse'};
invtype.values = {'GS', 'ARD', 'MSP', 'LOR', 'IID','EBB','EBBsparse'};
invtype.val = {'GS'};

woi = cfg_entry;
woi.tag = 'woi';
woi.name = 'Time window of interest';
woi.strtype = 'r';
woi.num = [1 2];
woi.val = {[-Inf Inf]};
woi.help = {'Time window to include in the inversion (ms)'};


foi = cfg_entry;
foi.tag = 'foi';
foi.name = 'Frequency window of interest';
foi.strtype = 'r';
foi.num = [1 2];
foi.val = {[0 256]};
foi.help = {'Frequency window (the same as high-pass and low-pass in the GUI)'};

hanning = cfg_menu;
hanning.tag = 'hanning';
hanning.name = 'PST Hanning window';
hanning.help = {'Multiply the time series by a Hanning taper to emphasize the central part of the response.'};
hanning.labels = {'yes', 'no'};
hanning.values = {1, 0};
hanning.val = {1};

patchfwhm = cfg_entry;
patchfwhm.tag = 'patchfwhm';
patchfwhm.name = 'Patch smoothness';
patchfwhm.strtype = 'r';
patchfwhm.num = [1 1];
patchfwhm.val = {[0.6]};
patchfwhm.help = {'Width of priors in cortex arb units (see inverse.smoothmm to see FWHM mm'};

npatches = cfg_entry;
npatches.tag = 'npatches';
npatches.name = 'Number of randomly selected patches';
npatches.strtype = 'i';
npatches.num = [1 1];
npatches.val = {[512]};
npatches.help = {'Number of randomly centred patches (priors) on each iteration'};

niter = cfg_entry;
niter.tag = 'niter';
niter.name = 'Number of iterations';
niter.strtype = 'i';
niter.num = [1 1];
niter.val = {[8]};
niter.help = {'Number of iterations'};


randpatch = cfg_branch;
randpatch.tag = 'randpatch';
randpatch.name = 'Random Patches';
randpatch.help = {'Define random patches'};
randpatch.val  = {npatches,niter};


% fixedpatch = cfg_entry;
% fixedpatch.tag = 'fixedpatch';
% fixedpatch.name = 'File with vertex indices for patch centres';
% fixedpatch.strtype = 's';
% %fixedpatch.num = [1 1];
% fixedpatch.val = {''};
% fixedpatch.help = {'Mat file with array Ip containing rows of iterations and columns of indices (patch centres for each iteration)'};

fixedfile = cfg_files;
fixedfile.tag = 'fixedfile';
fixedfile.name = 'Patch definition file ';
fixedfile.filter = 'mat';
fixedfile.num = [1 1];
fixedfile.help = {'Select patch definition file (mat file with variable Ip: rows are iterations columns are patch indices '};

fixedrows = cfg_entry;
fixedrows.tag = 'fixedrows';
fixedrows.name = 'Rows from fixed patch file';
fixedrows.strtype = 'i';
fixedrows.num = [1 2];
fixedrows.val = {[1 Inf]};
fixedrows.help = {'Select first and last row of patch sets to use from file'};


fixedpatch = cfg_branch;
fixedpatch.tag = 'fixedpatch';
fixedpatch.name = 'Fixed Patches';
fixedpatch.help = {'Define fixed patches'};
fixedpatch.val  = {fixedfile,fixedrows};



isfixedpatch = cfg_choice;
isfixedpatch.tag = 'isfixedpatch';
isfixedpatch.name = 'Patch definition';
isfixedpatch.help = {'Choose whether to use random or fixed patch centres.'};
isfixedpatch.values = {randpatch, fixedpatch};
isfixedpatch.val = {randpatch};


mselect = cfg_menu;
mselect.tag = 'mselect';
mselect.name = 'Selction of winning model';
mselect.help = {'How to get the final current density estimate from multiple iterations'};
mselect.labels = {'MaxF','Merge'};
mselect.values = {0,1};
mselect.val = {0};

nsmodes = cfg_entry;
nsmodes.tag = 'nsmodes';
nsmodes.name = 'Number of spatial modes';
nsmodes.strtype = 'i';
nsmodes.num = [1 1];
nsmodes.val = {[100]};
nsmodes.help = {'Number of spatial modes'};

umodes = cfg_files;
umodes.tag = 'umodes';
umodes.name = 'File containing spatial mode definition';
umodes.filter = 'mat';
umodes.num=[1 1];
umodes.val={''};
umodes.help = {'File (SPMU..) containing spatail mode definition, if no file is supplied bases these modes on lead fields'};


ntmodes = cfg_entry;
ntmodes.tag = 'ntmodes';
ntmodes.name = 'Number of temporal modes';
ntmodes.strtype = 'i';
ntmodes.num = [1 1];
ntmodes.val = {[4]};
ntmodes.help = {'Number of temporal modes'};


priorsmask  = cfg_files;
priorsmask.tag = 'priorsmask';
priorsmask.name = 'Priors file';
priorsmask.filter = '(.*\.gii$)|(.*\.mat$)|(.*\.nii(,\d+)?$)|(.*\.img(,\d+)?$)';
priorsmask.num = [0 1];
priorsmask.help = {'Select a mask or a mat file with priors.'};
priorsmask.val = {{''}};

space = cfg_menu;
space.tag = 'space';
space.name = 'Prior image space';
space.help = {'Space of the mask image.'};
space.labels = {'MNI', 'Native'};
space.values = {1, 0};
space.val = {1};

priors = cfg_branch;
priors.tag = 'priors';
priors.name = 'Source priors';
priors.help = {'Restrict solutions to pre-specified VOIs'};
priors.val  = {priorsmask, space};

locs  = cfg_entry;
locs.tag = 'locs';
locs.name = 'Source locations';
locs.strtype = 'r';
locs.num = [Inf 3];
locs.help = {'Input source locations as n x 3 matrix'};
locs.val = {zeros(0, 3)};

radius = cfg_entry;
radius.tag = 'radius';
radius.name = 'Radius of VOI (mm)';
radius.strtype = 'r';
radius.num = [1 1];
radius.val = {32};

restrict = cfg_branch;
restrict.tag = 'restrict';
restrict.name = 'Restrict solutions';
restrict.help = {'Restrict solutions to pre-specified VOIs'};
restrict.val  = {locs, radius};

custom = cfg_branch;
custom.tag = 'custom';
custom.name = 'Custom';
custom.help = {'Define custom settings for the inversion'};
custom.val  = {invfunc,invtype, woi, foi, hanning,isfixedpatch,patchfwhm,mselect,nsmodes,umodes,ntmodes, priors, restrict,outinv};

isstandard = cfg_choice;
isstandard.tag = 'isstandard';
isstandard.name = 'Inversion parameters';
isstandard.help = {'Choose whether to use standard or custom inversion parameters.'};
isstandard.values = {standard, custom};
isstandard.val = {standard};

crossval = cfg_entry;
crossval.tag = 'crossval';
crossval.name = 'Cross validation parameters';
crossval.strtype = 'r';
crossval.num = [1 2];
crossval.val = {[0 1]};
crossval.help = {'Percentage of cross validation test trials (eg 10 for 10 percent test trials, 0 for no cross val), number of cross val permutations. If number permutations*percentage equals 100, then make unique, non-repeating, selection of channels.'};

modality = cfg_menu;
modality.tag = 'modality';
modality.name = 'Select modalities';
modality.help = {'Select modalities for the inversion (only relevant for multimodal datasets).'};
modality.labels = {'All', 'EEG', 'MEG', 'MEGPLANAR', 'EEG+MEG', 'MEG+MEGPLANAR', 'EEG+MEGPLANAR'};
modality.values = {
    {'All'}
    {'EEG'}
    {'MEG'}
    {'MEGPLANAR'}
    {'EEG', 'MEG'}
    {'MEG', 'MEGPLANAR'}
    {'EEG', 'MEGPLANAR'}
    }';
modality.val = {{'All'}};

invert = cfg_exbranch;
invert.tag = 'invertiter';
invert.name = 'Source inversion, iterative';
invert.val = {D, val, whatconditions, isstandard, modality,crossval};
invert.help = {'Run imaging source reconstruction'};
invert.prog = @run_inversion;
invert.vout = @vout_inversion;
invert.modality = {'EEG'};

function  out = run_inversion(job)


D = spm_eeg_load(job.D{1});




inverse = [];
try
    inverse.priors=D.inv{job.val}.inverse.priors;
catch
    inverse.priors=[];
end;


if isfield(job.whatconditions, 'condlabel')
    inverse.trials = job.whatconditions.condlabel;
end
if numel(job.D)>1,
    error('iterative routine only meant for single subjects');
end;

savecopyinv=0; %% by default update dataset
if isfield(job.isstandard, 'custom')
    funccall=job.isstandard.custom.invfunc;
    inverse.type = job.isstandard.custom.invtype;
    inverse.woi  = fix([max(min(job.isstandard.custom.woi), 1000*D.time(1)) min(max(job.isstandard.custom.woi), 1000*D.time(end))]);
    inverse.Han  = job.isstandard.custom.hanning;
    inverse.lpf  =  fix(min(job.isstandard.custom.foi)); %% hpf and lpf are the wrong way round at the moment but leave for now
    inverse.hpf  =  fix(max(job.isstandard.custom.foi));
    inverse.mergeflag=job.isstandard.custom.mselect;
    
    savecopyinv = ~isempty(job.isstandard.custom.outinv);
    
    if ~isfield(job.isstandard.custom.isfixedpatch,'fixedpatch'), % fixed or random patch
        inverse.Np =  fix(max(job.isstandard.custom.isfixedpatch.randpatch.npatches));
        Npatchiter =  fix(max(job.isstandard.custom.isfixedpatch.randpatch.niter));
        allIp=[];
    else
        
        %% load in patch file
        
        dum=load(char(job.isstandard.custom.isfixedpatch.fixedpatch.fixedfile));
        if ~isfield(dum,'Ip'),
            error('Need to have patch indices in structure Ip');
        end;
        allIp=dum.Ip;
        lastrow=min(job.isstandard.custom.isfixedpatch.fixedpatch.fixedrows(2),size(allIp,1));
        rowind=[job.isstandard.custom.isfixedpatch.fixedpatch.fixedrows(1):lastrow];
        
        allIp=allIp(rowind,:);
        inverse.Np =  size(allIp,2);
        Npatchiter =  size(allIp,1);
        disp(sprintf('Using %d iterations of %d patches',Npatchiter,inverse.Np));
    end;
    inverse.Nm =  fix(max(job.isstandard.custom.nsmodes));
    
    inverse.Nt =  fix(max(job.isstandard.custom.ntmodes));
    if inverse.Nt==0,
        inverse.Nt=[];
        disp('Auto selecting number of modes');
    end;
    inverse.smooth=job.isstandard.custom.patchfwhm;
    testchans=[];
    crossU={};
    a=[];
    if ~isempty(job.isstandard.custom.umodes),
        disp('Loading spatial modes from file');
        
        a=load(cell2mat(job.isstandard.custom.umodes));
        
        if isfield(a,'testchans'), %% see if some of channels have been turned off
            testchans=a.testchans;
            if numel(a.U)~=size(testchans,1),
                error('U should have same number of elements as testchans has rows')
            end;
            crossU=a.U;
        else
            crossU={a.U};
        end;
    else %% if isempty U mode file
        A=[];
    end;
    
    
    
    
    
    if inverse.Nt==0,
        disp('Getting number of temporal modes from data');
        inverse.Nt=[];
    end;
    
    
    
    P = char(job.isstandard.custom.priors.priorsmask);
    if ~isempty(P)
        [p,f,e] = fileparts(P);
        switch lower(e)
            case '.gii'
                g = gifti(P);
                inverse.pQ = cell(1,size(g.cdata,2));
                for i=1:size(g.cdata,2)
                    inverse.pQ{i} = double(g.cdata(:,i));
                end
            case '.mat'
                load(P);
                inverse.pQ = pQ;
            case {'.img', '.nii'}
                S.D = D;
                S.fmri = P;
                S.space = job.isstandard.custom.priors.space;
                D = spm_eeg_inv_fmripriors(S);
                inverse.fmri = D.inv{D.val}.inverse.fmri;
                load(inverse.fmri.priors);
                inverse.pQ = pQ;
            otherwise
                error('Unknown file type.');
        end
    end
    
    if ~isempty(job.isstandard.custom.restrict.locs)
        inverse.xyz = job.isstandard.custom.restrict.locs;
        inverse.rad = job.isstandard.custom.restrict.radius;
    end
else %% standard inversion option, empty fields so they will become defaults
    funccall='Current';
    inverse.Np = [];
    inverse.Nt = [];
    inverse.Niter=[];
    inverse.Ns=[];
    Npatchiter = 1;
end

[mod, list] = modality(D, 1, 1);
if strcmp(job.modality{1}, 'All')
    inverse.modality  = list;
else
    inverse.modality  = intersect(list, job.modality);
end

if numel(inverse.modality) == 1
    inverse.modality = inverse.modality{1};
end

D = {};



for i = 1:numel(job.D)
    D{i} = spm_eeg_load(job.D{i});
    
    D{i}.val = job.val;
    
    D{i}.con = 1;
    
    if ~isfield(D{i}, 'inv')
        error('Forward model is missing for subject %d.', i);
    elseif  numel(D{i}.inv)<D{i}.val || ~isfield(D{i}.inv{D{i}.val}, 'forward')
        if D{i}.val>1 && isfield(D{i}.inv{D{i}.val-1}, 'forward')
            D{i}.inv{D{i}.val} = D{i}.inv{D{i}.val-1};
            warning('Duplicating the last forward model for subject %d.', i);
        else
            error('Forward model is missing for subject %d.', i);
        end
    end
    val=D{i}.val;
    D{i}.inv{val}.inverse = inverse;
    
    D{i}.inv{val}.inverse.allF=zeros(1,Npatchiter);
    
    
    %% commented out section will add an inversion at new indices
    
    %     for iterval=1:Npatchiter-1,
    %         D{i}.inv{iterval+val}=D{i}.inv{val}; %% copy inverse to all iterations which will be stored in the same file in higher vals
    %     end;
    D{i}.inv{D{i}.val}.inverse.PostMax=zeros(  length(D{i}.inv{D{i}.val}.forward.mesh.vert),1);
    
    
    M1.vertices  = D{i}.inv{val}.mesh.tess_mni.vert;
    M1.faces  = D{i}.inv{val}.mesh.tess_mni.face;
    
    
    
    
end


%% cross val set up




modalities = D{1}.inv{1}.forward.modality;
if size(modalities,1)>1,
    error('only works for single modality');
end;

megind=D{1}.indchantype(modalities);
origbadind=D{1}.badchannels;
megind=setdiff(megind,origbadind); %% start with just the good meg channels

Ntestchans=round(job.crossval(1)*length(megind)/100);
Nblocks=round(job.crossval(2));


if isfield(a,'testchans'),
    fprintf('Using testchans defined in spatial modes file\n')
    if Ntestchans~=size(testchans,2),
        error('Ntestchans different in U spec and cross val');
    end;
    if Nblocks>size(testchans,1),
        error('Not enough elements defined in U for this many iterations'),
    end;
else,
    fprintf('Making random test channels \n')
    rng('default');rng(0); %% always use same random seed here- so that cross vals (and Fs) between models can be compared
    
    for b=1:Nblocks,
        dum=randperm(length(megind));
        testchans(b,:)=dum(1:Ntestchans); %% channels which will be used for testing (Set to bad during training)
        %useind(b,:)=setdiff(1:length(megind),testchans(b,:)); %% training channels used here to get the spatial modes
    end;
    
end;





%% load in lead field matrix (with all good chans and make a copy)

[Lfull D{1}]=spm_eeg_lgainmat(D{1});
 gainname = D{1}.inv{1}.gainmat
% load(fullfile(D{1}.path, fname),'G','label');
% fullgainname=fullfile(D{1}.path, ['Full_' fname]);
% save(fullgainname,'G','label');

Nchans=D{1}.nchannels;

crosserr=zeros(Nblocks,1);
allrms=zeros(Nblocks,1);
crossF=zeros(Nblocks,1);
xvalchans=testchans.*0;
for b=1:Nblocks,
    
    
    badmegind=testchans(b,:);
    xvalchans(b,:)=megind(badmegind);
    D{1} = badchannels(D{1}, 1:Nchans, 0); %% set all chans to good
    D{1} = badchannels(D{1}, [megind(badmegind) origbadind], 1); %% set orignal + test chans as bad
    
    
    
    D{1}.inv{val}.inverse.A=crossU{b};;
    
    [Dout,allmodels,allF] = spm_eeg_invertiter(D,Npatchiter,funccall,allIp);
    
    Dout{1} = badchannels(Dout{1}, 1:Nchans, 0); %% set all channels as good
    Dout{1} = badchannels(Dout{1}, origbadind, 1); %% set orginal bad channels back to bad
    inverse=Dout{1}.inv{val}.inverse;
    
    
    
    L=inverse.L; %% lead fields used
    U=inverse.U{:}; %% lead field projection
    T=inverse.T; %% temporal projector
    
    Ic=inverse.Ic{1}; % trials
    It=inverse.It; %% time points
    Ik=inverse.Ik; %% trials
    
    M=inverse.M; %% MAP operator
    
    
    crossF(b)=inverse.F; %% ok now trained for this model
    
    %% use hanning window if used in inversion
    if inverse.Han,
        W  = repmat(spm_hanning(length(It))',length(Ic),1); %% use hanning unless specified
    else
        W=ones(length(It),length(Ic));
    end;
    
    
    if Ntestchans>0, %% if there are any channels to predict..
        fprintf('\n Cross val iteration %d of %d, first two chans:%d, %d..\n',b,Nblocks,testchans(b,1),testchans(b,2));
        errpred=0;
        rmstot=0;
        
        fprintf('Running cross val..');
        for t=1:length(Ik),
            traindata=squeeze((D{1}(Ic,It,Ik(t))));
            
            traindata=traindata.*W; %% hanning window
            allUYtrain=U*(traindata)*T; %
            Jtrain=M*allUYtrain; %% source estimate based on training data
            Ypred=Lfull*Jtrain; %% now predict test channels based on trained source estimate
            Ymeas=squeeze((D{1}(megind(testchans(b,:)),It,Ik(t)))); % test channels - not used in inversion
            Ymeas=Ymeas.*W(1:size(testchans,2),:)*T;
            errpred=errpred+sqrt(mean((sum(Ypred(badmegind,:)-Ymeas).^2))); % asusming same model in trial and test (i.e. not necessarily time locked)
            rmstot=rmstot+sqrt(mean(sum(Ymeas.^2)));
            
        end;
        crosserr(b)=errpred;
        allrms(b)=rmstot;
    end; %if pctest>0
    
    
end; % for b



D{1}.inv{val}.inverse=inverse;
D{1}.inv{val}.inverse.crosserr=crosserr;
D{1}.inv{val}.inverse.crossF=crossF;
D{1}.inv{val}.inverse.allrms=allrms;
D{1}.inv{val}.inverse.xvalchans=xvalchans;

D{1} = badchannels(D{1}, 1:Nchans, 0); %% set all chans to good
D{1} = badchannels(D{1}, [origbadind], 1); %% set orignal bad chans as bad

if ~iscell(D)
    D = {D};
end

for i = 1:numel(D)
    save(D{i});
    if savecopyinv, %% actually going to save the dataset anyway- badchannels need to be set right
        outinvname=[D{i}.path filesep job.isstandard.custom.outinv '_' D{i}.fname];
        fprintf(sprintf('\nNB. Original dataset is being updated and a copy of inversion results written to %s',outinvname));
        
        inv=D{i}.inv{val};
        spmfilename=D{i}.fname;
        save(outinvname,'-v7.3','inv','spmfilename');
    end;
end

out.D = job.D;

function dep = vout_inversion(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'M/EEG dataset(s) after imaging source reconstruction';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});

