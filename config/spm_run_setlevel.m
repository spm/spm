function out = spm_run_setlevel(job)
% out = spm_run_setlevel(job)
% test to see how likely it is that an SPM statistical image is a random field.
% based on:
%  Set-level threshold-free tests on the intrinsic volumes of SPMs.
%   Barnes GR, Ridgway GR, Flandin G, Woolrich M, Friston K. Neuroimage. 2013
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% $Id: spm_run_setlevel.m 5554 2013-06-13 09:13:56Z gareth $

%-Load SPM.mat file
%--------------------------------------------------------------------------
max_resel_dimension=3;  %% for now
LOWTHRESH=7; %% t value
threshlevels=[-LOWTHRESH:0.2:LOWTHRESH];
PLOTON=0;
ESTZEROLKC=0; %% don't try and estimate the 0th LKC (as this is given by the topology of the mask)

SPM = [];
load(job.spmmat{:});
out.spmmat = job.spmmat;

cindex=job.cindex; %% contrast of interest

%% First we need to redo the residual images - these have normally been cleaned up.

Ic=NaN; %% adjust for everything
disp('Writing new residual images');
Vres = spm_write_residuals(SPM,Ic);

%%% Now to make standarized residual images
nSres=size(Vres,2); %% number of subjects/observations
DIM=Vres.dim;
%-Initialise standardised residual images
%--------------------------------------------------------------------------

M=Vres.mat;
VResI(1:nSres) = deal(struct(...
    'fname',   [],...
    'dim',     DIM,...
    'dt',      [spm_type('float64') spm_platform('bigend')],...
    'mat',     M,...
    'pinfo',   [1 0 0]',...
    'descrip', 'spm_spm:StandardisedResiduals'));

for i = 1:nSres
    VResI(i).fname   = [sprintf('ResI_%04d', i) spm_file_ext];
    VResI(i).descrip = sprintf('spm_spm:ResI (%04d)', i);
end
VResI = spm_data_hdr_write(VResI);


ResMS=spm_read_vols(SPM.VResMS); %% residual mean square ( ResSS divided by tr(RV))

ResSS=ResMS.*SPM.xX.trRV; %%  residual sum of squares
%% now write out standardized residuals for each subject/observation
outvol            = NaN(size(ResSS));
cmask=find(isfinite(ResSS));

for j=1:nSres
    res=spm_read_vols(Vres(j));
    outvol(cmask) = res(cmask)./(sqrt(ResSS(cmask)/SPM.xX.erdf)); % or xX.erdf
    VResI(j) = spm_data_write(VResI(j), outvol);
end

maskV=SPM.VM;

mask=spm_read_vols(maskV);
maskind=intersect(cmask,find(mask));


datatrial1d=zeros(nSres,length(maskind));

%% disp read back in standardised residuals
for t=1:nSres,
    [sigvoldata,XYZ]=spm_read_vols(VResI(t));
    datatrial1d(t,:)=sigvoldata(maskind);
end % for t


%% read in statistical image
[Teststat,XYZ]=spm_read_vols(SPM.xCon(cindex).Vspm);
test_STAT=SPM.xCon(cindex).STAT;


ec_R0=zeros(nSres+1,1); %% topology
alleuler2_spm=zeros(nSres+1,length(threshlevels));

disp(sprintf('Getting ECs over %d residual images and %d thresholds...',nSres,length(threshlevels)));

for tr=1:nSres+1, % last extra trial is for Teststat
    vol_data=ones(size(sigvoldata)).*NaN;
    if tr<=nSres,
        vol_data(maskind)=datatrial1d(tr,:);
    else
        vol_data(maskind)=Teststat(maskind);
    end;
    
    binpoints=zeros(size(vol_data));
    binpoints(find(vol_data))=1;
    R0 = spm_resels_vol(binpoints,SPM.xVol.FWHM);
    ec_R0(tr)=R0(1);
    %% now get EC at each threshold and each observation
    for threshind=1:length(threshlevels),
        
        useind=find(vol_data>=threshlevels(threshind));
        binpoints=zeros(size(vol_data));
        binpoints(useind)=1;
        R0 = spm_resels_vol(binpoints,SPM.xVol.FWHM);
        ec_spm=R0(1);
        ECcount(tr,threshind)=ec_spm;
    end; % for threshind
end; % for tr

if max(ec_R0)~=min(ec_R0)
    error('subjects have different topology');
end;
R0_set=ec_R0(1);


%% GET EC DENSITIES FOR INDIVIDUAL TRIAL RESIDUALS (Z DIST) AND FOR TEST OF INTEREST (independent of data)
allpju_trial=[]; %% independent of data
allpju_test=[]; %% independent of data
df=[1 nSres-size(SPM.xX.X,2)];
trial_STAT='Z';
fitind=find(abs(threshlevels)<=LOWTHRESH);

for threshind=1:length(fitind),
    %%% NEED TO CHECK THIS ___
    [ECperresel_trial]=spm_ECdensity(trial_STAT,threshlevels(fitind(threshind)),df); %% ACTUALLY THIS IS LKC density
    [ECperresel_test]=spm_ECdensity(test_STAT,threshlevels(fitind(threshind)),df);
    
    for d=0:3, %% Turn density estimates from EC per resel to EC per LKC unit
        ECperLKC_trial(d+1,:)=ECperresel_trial(d+1,:)./power(4*log(2),d/2);
        ECperLKC_test(d+1,:)=ECperresel_test(d+1,:)./power(4*log(2),d/2);
    end;
    
    allpju_trial(threshind,:)=ECperLKC_trial(1:max_resel_dimension+1)';
    allpju_test(threshind,:)=ECperLKC_test(1:max_resel_dimension+1)';
end;

%%% GET RESEL ESTIMATES BASED ON SMOOTHNESS
[RPV]=spm_read_vols(SPM.xVol.VRpv); %% resel per voxel image
allR=SPM.xVol.R; %% resel counts





reselVec=SPM.xVol.R;
LKCresel=[];
for d=0:max_resel_dimension, %%% or alternatively LKC estimate from the reselts
    LKCresel(d+1)=reselVec(d+1)*power(4*log(2),d/2);
end;



for tbase=1:nSres+2,
    allpju=[];
    Y=[];
    if tbase<=nSres,
        epochtype=0; %% residual from trial
    else
        if tbase==nSres+1, %% test statistic
            epochtype=1;
        else
            epochtype=2; %% mean of ECs over trials
        end;
    end;
    
    
    
    
    switch epochtype,
        case 0, %% single residual trial
            Y=squeeze(ECcount(tbase,fitind))'; %% get LKC based on single residual image for these N trials
            allpju=allpju_trial; %% density for trial
        case 1 % single t stat
            Y=squeeze(ECcount(nSres+1,fitind))'; %% get LKC based on single Teststat image for these N trials
            allpju=allpju_test; %% density for test
        case 2, %% mean of all trials
            Y=squeeze(mean(ECcount(1:nSres,fitind),1))'; %% get LKC based on average EC over trials in iteration k
            %allaverageY(k,:)=Y;
            allpju=allpju_trial; %% density for trial
    end;
    
    
   
    dimension_test=max_resel_dimension;
    
    %% LKC based on average EC through basic regression
    if ~ESTZEROLKC, %% do not estimate 0th LKC
        LKC0=R0_set; %% TAKE THIS AS FIXED
        Ydash=Y-LKC0*allpju(:,1);
        useLKCind=2:dimension_test+1; 
        LKC_est=pinv(allpju(:,useLKCind))*Ydash; %% get least squares estimate of extra LKC coeffs
        LKC=[LKC0; LKC_est];  %% put them back together with 0th order
       
    else %% estimate allLKCs
        useLKCind=1:dimension_test+1;
        LKC=pinv(allpju(:,useLKCind))*Ydash; %% get least squares estimate of all LKC coeffs
    end; % if
    
    
    
    allLKCregress(tbase,:)=LKC; %% LKC based on fit to average EC (/ or EC of Teststat)
    
end; % for tbase



gY = [squeeze(allLKCregress(nSres+1,useLKCind)); squeeze(allLKCregress(1:nSres,useLKCind))];
%then the design matrix and contrast could be e.g.
gX = [1 0; [zeros(nSres, 1) ones(nSres, 1)]];
gC = [1 -1]';
[CVA] = spm_cva(gY,gX,[],gC); %% do multivariate test


%% get empirical mean and sd of EC over subjects / observations.

meanECtest_regtrial=mean(squeeze(allLKCregress(1:nSres,:)))*allpju_test'; %% unweighted based on mean Euler
sdECtest_regtrial=std(squeeze(allLKCregress(1:nSres,:)))*allpju_test'; %% unweighted based on mean Euler

%% estimate what EC profile should be based on smoothness of image
meanECtest_resel=LKCresel*allpju_test';

%% PLOT RESULTS

Fgraph  = spm_figure('GetWin','Graphics'); spm_figure('Clear',Fgraph);

h=plot(threshlevels,squeeze(ECcount(nSres+1,:,:)),'r-',threshlevels,meanECtest_regtrial,'b:',threshlevels,meanECtest_resel,'go');
set(h,'LineWidth',3);
hold on;
errorbar(threshlevels,meanECtest_regtrial,sdECtest_regtrial,'.b');
set(gca,'Fontsize',18);
set(gcf,'color','w');
set(h,'LineWidth',3);
hold on;

xlabel('threshold');
ylabel('EC');
legend(sprintf('Observered EC for %s field',test_STAT),...
    sprintf('Random %s field based on regression',test_STAT),...
    sprintf('Random %s field based on smoothness',test_STAT));



title(sprintf('Probability that this is a random field  p<%3.4f ',CVA.p));



%%% FINISHED

%-Move to the directory where the SPM.mat file is
%--------------------------------------------------------------------------
original_dir = pwd;
cd(fileparts(job.spmmat{:}));
cd(original_dir);

fprintf('Done\n')
return
