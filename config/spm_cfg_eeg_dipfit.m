function dipfit = spm_cfg_eeg_dipfit
% Configuration file for configuring imaging source inversion
% reconstruction
%__________________________________________________________________________
% Copyright (C) 2010-2016 Wellcome Trust Centre for Neuroimaging

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
locs.name = 'Prior Source locations';
locs.strtype = 'r';
locs.num = [Inf 3];
locs.help = {'Input source locations as n x 3 matrix'};
locs.val = {zeros(0, 3)};

locvar  = cfg_entry;
locvar.tag = 'locvar';
locvar.name = 'Prior Source location variance';
locvar.strtype = 'r';
locvar.num = [Inf 3];
locvar.help = {'Source location variance as n x 3 matrix. So [100 100 100] means approx 10mm uncertainty in any direction'};
locvar.val = {ones(1, 3)*100};

moms  = cfg_entry;
moms.tag = 'moms';
moms.name = 'Prior Source moments';
moms.strtype = 'r';
moms.num = [Inf 3];
moms.help = {'Input source moments as n x 3 matrix. If unsure of orientation leave this as zeros'};
moms.val = {zeros(1, 3)};

momvar  = cfg_entry;
momvar.tag = 'momvar';
momvar.name = 'Prior Source moment variance';
momvar.strtype = 'r';
momvar.num = [Inf 3];
momvar.help = {'Source moment variance as n x 3 matrix. So [100 100 100] means that expect to have magnitudes of around 10nAm'};
momvar.val = {ones(1, 3)*100};

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

dipfit = cfg_exbranch;
dipfit.tag = 'dipfit';
dipfit.name = 'Bayesian Dipole fit';
dipfit.val = {D, val, whatconditions, woi, locs, locvar,moms,momvar, niter,modality};
dipfit.help = {'Run imaging source reconstruction'};
dipfit.prog = @run_dipfit;
dipfit.vout = @vout_dipfit;
%dipfit.modality = {'MEG'};

function  out = run_dipfit(job)




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
    
end
D=D{1};
megind=D.indchantype(job.modality);
megind=setdiff(megind,D.badchannels);

usesamples=intersect(find(D.time>=job.woi(1)./1000),find(D.time<=job.woi(2)./1000));

    
condind=strmatch(job.whatconditions.condlabel,D.conditions)


% if length(condstr)>1,
%     error('not implemented for more than 1 condition yet');
% end;

fitdata=D(megind,usesamples,condind);


pos=coor2D(D)';                                                        %Uses info in MEG datafile (D) to print list of 2D locations of each channel (pos)
labels=num2str(megind');                                               %Prints numerical label of each channel (labels)




s0_mni=job.locs;
w0_mni=job.moms;
diags_s0_mni=job.locvar;
diags_w0_mni=job.momvar;

ndips=size(job.locs,1); %% single dipole
%% Set up the forward model
val=D.val;                                                             %Use the most recent forward model saved in D
P=[];
P.y=mean(fitdata,2);                                                    %Z=avdata at each sensor at 'fitind' (timepoint identified in VB_ECD_FindPeak)
P.forward.sens = D.inv{val}.datareg.sensors;
P.forward.vol  = D.inv{val}.forward.vol;
P.modality = D.modality;
P.Ic = megind;
P.channels = D.chanlabels(P.Ic);

spm_eeg_plotScalpData(P.y,pos,labels);
%% Transform to mni space
M1 = D.inv{val}.datareg.toMNI;
[U, L, V] = svd(M1(1:3, 1:3));
orM1(1:3,1:3) =U*V';                                                   %For switching orientation between meg and mni space


mnivol = ft_transform_vol(M1, P.forward.vol);                          %Used for inside head calculation


%% Put priors in ctf space
alldiag_Wcov=[];
alldiag_Scov=[];

figure;
hold on;

for d=1:ndips,
    
    ind=(d-1)*3+1:d*3;
    Priors.mu_w0(ind)= orM1*w0_mni(ind)';
    
    ctfpos=D.inv{val}.datareg.fromMNI*[s0_mni(ind) 1]';
    plot3(ctfpos(1)./1000,ctfpos(2)./1000,ctfpos(3)./1000,'r*');
    Priors.mu_s0(ind)=ctfpos(1:3); % ./1000;
    alldiag_Wcov(ind)=ones(3,1)*diags_w0_mni(d);
    alldiag_Scov(ind)=ones(3,1)*diags_s0_mni(d);
end;



Priors.mu_s0=spm_vec(Priors.mu_s0);
Priors.mu_w0=spm_vec(Priors.mu_w0);
Priors.S_s0=diag(alldiag_Scov);
Priors.S_w0= diag(alldiag_Wcov);


P.priors=Priors;

%% Calculate the covariance of the raw data

%% Next line calculates the SNR at each sensor, for the peak you've identified

%% The SNRamp is the average SNR for the peak you've identified across all sensors
SNRamp=3;

%% So, the SNR obviously varies across the brain, and is different at each sensor. How much of the variance observed in the sensor data is due to noise, and how much is due to real signal?
%P.priors.hE=log(SNRamp^2);                                             %Expected log precision of data
P.priors.hE=log(SNRamp);  %% seemed to be closets to value arrived at in line serach                                            %Expected log precision of data
%% run a line search over hE values
%hevals=[0.001 0.1 0.5 1 2] %  1 2 5 10]

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
    P.forward.vol = ft_read_vol(P.forward.vol);
end
P.forward.sens    = data.(P.modality(1:3)).sens;
P.forward.siunits = data.siunits;
% M1 = data.transforms.toMNI;
% [U, L, V] = svd(M1(1:3, 1:3));
% orM1(1:3,1:3) = U*V';                                                  %For switching orientation between meg and mni space

if isempty(P.Ic)
    error(['The specified modality (' P.modality ') is missing from file ' D.fname]);
else
    P.channels = D.chanlabels(P.Ic);
end
P.forward.chanunits = D.units(P.Ic);
[P.forward.vol, P.forward.sens] =  ft_prepare_vol_sens( ...
    P.forward.vol, P.forward.sens, 'channel', P.channels);
plot3(P.forward.sens.chanpos(:,1),P.forward.sens.chanpos(:,2),P.forward.sens.chanpos(:,3),'bo');

for j=1:job.niter;
    Pout(j)        = spm_eeg_inv_vbecd(P);
    varresids(j)   = var(Pout(j).y-Pout(j).ypost);
    pov(j)         = 100*(1-varresids(j)/var(Pout(j).y));          %Percent variance explained
    allF(j)        = Pout(j).F;
    dip_mom        = reshape(Pout(j).post_mu_w,3,length(Pout(j).post_mu_w)/3);
    dip_amp(j,:)   = sqrt(dot(dip_mom,dip_mom));
    megloc         = reshape(Pout(j).post_mu_s,3,length(Pout(j).post_mu_s)/3); %Loc of dip (3 x n_dip)
    mniloc{j}      = D.inv{val}.datareg.toMNI*[megloc;ones(1,size(megloc,2))]; %Actual MNI location (with scaling)
    megmom{j}      = reshape(Pout(j).post_mu_w,3,length(Pout(j).post_mu_w)/3); %Moments of dip (3 x n_dip)
end                                                                %For j in serial loop

%%%% Save VB-ECD results into Output structure and overwrite datafile

[maxF,maxind]=max(allF);


inverse=[];
inverse.F=maxF;
inverse.modelmniloc=mniloc{maxind}(1:3)';
inverse.modelmnimom=megmom{maxind}(1:3)';
inverse.Pout=Pout(maxind);
inverse.P=P;
D.inv{job.val}.inverse=inverse;

hf = spm_figure('FindWin','Graphics');
figure(hf);
clf;
subplot(3,2,1);
plot(inverse.Pout.y,inverse.Pout.ypost,'o',inverse.Pout.y,inverse.Pout.y,':');
xlabel('measured');ylabel('modelled');
title(sprintf('Free energy=%3.2f',inverse.F));



axesY = axes(...
    'Position',[0.1 0.4 0.3 0.2],...
    'hittest','off');
in.f = hf;
in.noButtons = 1;
in.ParentAxes = axesY;


spm_eeg_plotScalpData(P.y,D.coor2D',P.channels,in)
title(axesY,'measured data')

axesY = axes(...
    'Position',[0.5 0.4 0.3 0.2],...
    'hittest','off');

in.ParentAxes=axesY;
spm_eeg_plotScalpData(inverse.Pout.ypost,D.coor2D',P.channels,in)
title(axesY,'Modelled data');

subplot(3,3,7); hold on;

M=D.inv{val}.mesh.tess_mni;
h=trisurf(M.face,M.vert(:,1),M.vert(:,2),M.vert(:,3));
set(h,'Facecolor','cyan');
set(h,'Edgecolor','none');
alpha(0.1)
plot3(inverse.modelmniloc(:,1),inverse.modelmniloc(:,2),inverse.modelmniloc(:,3),'k*');
view([0 0 1]);


subplot(3,3,8);  hold on;
M=D.inv{val}.mesh.tess_mni;
h=trisurf(M.face,M.vert(:,1),M.vert(:,2),M.vert(:,3));
set(h,'Facecolor','cyan');
set(h,'Edgecolor','none');
alpha(0.1)
plot3(inverse.modelmniloc(:,1),inverse.modelmniloc(:,2),inverse.modelmniloc(:,3),'k*');

view([0 1 0]);

subplot(3,3,9); hold on;
M=D.inv{val}.mesh.tess_mni;
h=trisurf(M.face,M.vert(:,1),M.vert(:,2),M.vert(:,3));
set(h,'Facecolor','cyan');
set(h,'Edgecolor','none');
alpha(0.1)
plot3(inverse.modelmniloc(:,1),inverse.modelmniloc(:,2),inverse.modelmniloc(:,3),'k*');

view([1 0 0]);

subplot(3,3,7);
text(-50,-110,'MNI coord:');
text(-50,-130,num2str(round(inverse.modelmniloc)));

subplot(3,3,9);
text(0,-190,-40,'MNI mom:');
text(0,-140,-50,num2str(round(inverse.modelmnimom)));

if ~iscell(D)
    D = {D};
end


for i = 1:numel(D)
    save(D{i});
end

out.D = job.D;

function dep = vout_dipfit(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'M/EEG dataset(s) after imaging source reconstruction';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'filter','mat'}});

