function momentfit = spm_cfg_eeg_momentfit
% Configuration file for imaging source inversion reconstruction.
% This version to supply position and orientation parameters idea is to
% estimate dipole moments given priors and return a model evidence for
% these priors.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_cfg_eeg_momentfit.m 7763 2020-01-02 15:01:38Z gareth $


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

all = cfg_const;
all.tag = 'all';
all.name = 'All';
all.val  = {1};
all.help = {''};

condlabel = cfg_entry;
condlabel.tag = 'condlabel';
condlabel.name = 'Condition label';
condlabel.strtype = 's';
condlabel.val = {''};
condlabel.help = {''};

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
whatconditions.help = {'What conditions to include?'};

whichgenerator = cfg_menu;
whichgenerator.tag = 'whichgenerator';
whichgenerator.name = 'Fit current dipole (in volume conductor) or magnetic dipole (in free space) ';
whichgenerator.labels = {'current',  'magnetic'};
whichgenerator.values = {1,0};
whichgenerator.help = {'Fitting field due to an electrical current within the head or a coil (magnetic dipole) outside of the head ?'};

woi = cfg_entry;
woi.tag = 'woi';
woi.name = 'Time window of interest';
woi.strtype = 'r';
woi.num = [1 2];
woi.val = {[-Inf Inf]};
woi.help = {'Time window to include in the inversion (ms)'};


hanning = cfg_menu;
hanning.tag = 'hanning';
hanning.name = 'PST Hanning window';
hanning.help = {'Multiply the time series by a Hanning taper to emphasize the central part of the response.'};
hanning.labels = {'yes', 'no'};
hanning.values = {1, 0};
hanning.val = {1};


locs  = cfg_entry;
locs.tag = 'locs';
locs.name = 'Source locations';
locs.strtype = 'r';
locs.num = [Inf 3];
locs.help = {'Source locations as n x 3 matrix'};
locs.val = {zeros(0, 3)};

source_ori  = cfg_entry;
source_ori.tag = 'source_ori';
source_ori.name  = 'Source orientations as n x 3 matrix';
source_ori.strtype = 'r';
source_ori.num  = [Inf 3];
source_ori.help = {'Source orientation as n x 3 matrix of unit vectors'};
source_ori.val = {ones(1, 3)};

moms  = cfg_entry;
moms.tag = 'moms';
moms.name = 'Prior mean of source moments';
moms.strtype = 'r';
moms.num = [Inf 1];
moms.help = {'Input source moments as n x 1 matrix (in nAm)'};
moms.val = {zeros(1, 1)};

momvar  = cfg_entry;
momvar.tag = 'momvar';
momvar.name = 'Prior variance of source moments';
momvar.strtype = 'r';
momvar.num = [Inf 1];
momvar.help = {'Prior source moment variance as n x 1 matrix (in (nAm) ^2)'};
momvar.val = {ones(1,1 )*100};

ampsnr  = cfg_entry;
ampsnr.tag = 'ampsnr';
ampsnr.name = 'Estimated amplitude SNR';
ampsnr.strtype = 'r';
ampsnr.num = [ 1 1];
ampsnr.help = {'Estimate amplitude SNR as a ratio (1 for equal signal and noise amplitude)'};
ampsnr.val = {10};

niter = cfg_entry;
niter.tag = 'niter';
niter.name = 'Number of fit iterations';
niter.strtype = 'n';
niter.help = {'Number of restarts of algorithm to avoid local extrema'};
niter.val = {10};

modality = cfg_menu;
modality.tag = 'modality';
modality.name = 'Select modalities';
modality.help = {'Select modalities for the inversion (only relevant for multimodal datasets).'};
modality.labels = {'EEG', 'MEG', 'EEG+MEG'};
modality.values = {
    {'EEG'}
    {'MEG'}
    {'EEG', 'MEG'}
    }';
modality.val = {{'MEG'}};

momentfit = cfg_exbranch;
momentfit.tag = 'momentfit';
momentfit.name = 'Bayesian moment fit';
momentfit.val = {D, val,whichgenerator, whatconditions, woi, locs, source_ori,moms,momvar,ampsnr, niter,modality};
momentfit.help = {'Run imaging source reconstruction'};
momentfit.prog = @run_momentfit;
momentfit.vout = @vout_momentfit;
%dipfit.modality = {'MEG'};

function  out = run_momentfit(job)




D = {};


for i = 1:numel(job.D)
    D{i} = spm_eeg_load(job.D{i});
    
    D{i}.val = job.val;
    
    D{i}.con = 1;
    
    if job.whichgenerator==1, %% if current dipole check for a forward model
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
    else
        disp('Magdip fit: Not checking for volume conductor (assuming infinite medium)');
    end; %% if whichgenerator
    
end

D=D{1};
megind=D.indchantype(job.modality);
megind=setdiff(megind,D.badchannels);

usesamples=intersect(find(D.time>=job.woi(1)./1000),find(D.time<=job.woi(2)./1000));


if isfield(job.whatconditions,'all') %% all conditions
    condind=setdiff(1:D.ntrials,D.badtrials);
else
    condind=strmatch(job.whatconditions.condlabel,D.conditions)
end;


fitdata=zeros(length(megind),length(usesamples)); %% pool over conditions
for f=1:length(condind),
    fitdata=fitdata+squeeze(D(megind,usesamples,condind(f)));
end;
fitdata=fitdata./length(condind);


pos=coor2D(D)';                                                        %Uses info in MEG datafile (D) to print list of 2D locations of each channel (pos)
labels=num2str(megind');                                               %Prints numerical label of each channel (labels)




dippos_mni=job.locs;
dipor_mni=job.source_ori;
priormom=job.moms;
priormomvar=job.momvar;


ndips=size(dippos_mni,1); %% single dipoles
if size(dipor_mni,1)~=ndips || size(priormom,1)~=ndips || size(priormomvar,1)~=ndips,
    error('Inputs for prior mean and uncertainty should all have same number of rows (= number dipoles)');
end;
%% Set up the forward model
val=D.val;                                                             %Use the most recent forward model saved in D
P=[];
P.y=mean(fitdata,2);                                                    %Z=avdata at each sensor at 'fitind' (timepoint identified in VB_ECD_FindPeak)
P.forward.sens = D.inv{val}.datareg.sensors;

P.modality = D.modality;
P.Ic = megind;
P.channels = D.chanlabels(P.Ic);

spm_eeg_plotScalpData(P.y,pos,labels);
%% Transform from MNI space
M1 = D.inv{val}.datareg.fromMNI;
[U, L, V] = svd(M1(1:3, 1:3));

orM1(1:3,1:3) =U*V';                                                   %For switching orientation between meg and mni space


% Mnivert=[D.inv{val}.mesh.tess_mni.vert ones(size(D.inv{val}.mesh.tess_mni.vert,1),1)];
% meshctf.vertices=D.inv{val}.datareg.fromMNI*Mnivert';
% meshctf.vertices=meshctf.vertices(1:3,:)';
% meshctf.faces=D.inv{val}.mesh.tess_mni.face;
% meshctfnorms=spm_mesh_normals(meshctf);
% meshctfnorms([6157 2087],:)
% mnivol = ft_transform_vol(M1, P.forward.vol);                          %Used for inside head calculation

figure;
hold on;

dippos_ctf=dippos_mni.*0; %% coords in MEG space
dipor_ctf=dippos_ctf;
for d=1:ndips,
    
    dipor_mni(d,:)=dipor_mni(d,:)./sqrt(dot(dipor_mni(d,:),dipor_mni(d,:)));; % force to unit vector
    dipor_ctf(d,:)= orM1*dipor_mni(d,:)'; %% want orientation in ctf space
    
    pos_ctf=D.inv{val}.datareg.fromMNI*[dippos_mni(d,:) 1]';
    dippos_ctf(d,:)=pos_ctf(1:3);
end;





%% The SNRamp is the average SNR for the peak you've identified across all sensors
SNRamp=job.ampsnr;

%P.priors.hE=log(SNRamp^2);                                             %Expected log precision of data
P.priors.hE=log(SNRamp);
%% run a line search over hE values

%P.priors.hE=hevals(hind); %% log(SNRamp^2);                                             %Expected log precision of data
hcpval=-2; %; %% was -2
P.priors.hC=exp(hcpval);                                      %varSNR; % variability of the above precision

%% This is the range of possible SNR values at each sensor that we're willing to expect. If the Hcpval is small, the range will be small.
disp('approx amp SNR range:')
exp(P.priors.hE/2-2*sqrt(P.priors.hC))
exp(P.priors.hE/2+2*sqrt(P.priors.hC))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Start the VB-ECD iterations for the ith subject


data = spm_eeg_inv_get_vol_sens(D, val, [], 'inv', P.modality);
P.forward.vol     = data.(P.modality(1:3)).vol;
if ischar(P.forward.vol)
    P.forward.vol = ft_read_headmodel(P.forward.vol);
end
if job.whichgenerator==0, %% if magnetic dipole override forward model
    
    P.forward.vol.type='infinite_magneticdipole';
end;
P.forward.sens    = data.(P.modality(1:3)).sens;
P.forward.siunits = data.siunits;

if isempty(P.Ic)
    error(['The specified modality (' P.modality ') is missing from file ' D.fname]);
else
    P.channels = D.chanlabels(P.Ic);
end
P.forward.chanunits = D.units(P.Ic);
[P.forward.vol, P.forward.sens] =  ft_prepare_vol_sens( ...
    P.forward.vol, P.forward.sens, 'channel', P.channels);
plot3(P.forward.sens.chanpos(:,1),P.forward.sens.chanpos(:,2),P.forward.sens.chanpos(:,3),'bo');

P.dippos_ctf=dippos_ctf;
P.dipor_ctf=dipor_ctf;
P.priors.mom=priormom;
P.priors.momvar=priormomvar;
for j=1:job.niter,
    
    Pout(j)        = spm_eeg_inv_vbecd_mom(P);
    varresids(j)   = var(Pout(j).y-Pout(j).ypost);
    pov(j)         = 100*(1-varresids(j)/var(Pout(j).y));          %Percent variance explained
    allF(j)        = Pout(j).F;
    
    
end                                                                %For j in serial loop

%%%% Save VB-ECD results into Output structure and overwrite datafile

[maxF,maxind]=max(allF);


inverse=[];

inverse.F=maxF;
inverse.Pout=Pout(maxind);
inverse.P=P;
D.inv{job.val}.inverse=inverse;


hf=figure;
plot(inverse.Pout.y,inverse.Pout.ypost,'o',inverse.Pout.y,inverse.Pout.y,':');
xlabel('measured');ylabel('modelled');
title(sprintf('Free energy=%3.2f',inverse.F));






if ~iscell(D)
    D = {D};
end


for i = 1:numel(D)
    save(D{i});
end

out.D = job.D;

function dep = vout_momentfit(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'M/EEG dataset(s) after imaging source reconstruction';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});

