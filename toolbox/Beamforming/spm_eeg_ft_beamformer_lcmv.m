function [stats,mnipositions]=spm_eeg_ft_beamformer_lcmv(S)
% Compute power-based beamformer image
% FORMAT [stats,mnipositions]=spm_eeg_ft_beamformer_lcmv(S)
%
% returns a stats structure containing univariate t test on power (based
% purely on sign in first column of design matrix S.design.X)
% and a list of the image files produced
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_eeg_ft_beamformer_lcmv.m 4344 2011-06-07 14:58:43Z gareth $

[Finter,Fgraph] = spm('FnUIsetup','univariate LCMV beamformer for power', 0);
%%

%% ============ Load SPM EEG file and verify consistency
if nargin == 0
    S = [];
end


if ~isfield(S,'gridpos'),
    S.gridpos=[];
    end;

 if ~isfield(S,'maskgrid'),
     %%
    S.maskgrid=[];
    end;

 
if ~isfield(S,'design'),
   error('Design matrix required');
end; % if

if ~isfield(S,'return_weights')
    ctf_weights=[];
    S.return_weights=0;
end

if ~isfield(S,'Niter')
    S.Niter=[];
end; % if 

if isempty(S.Niter),
    S.Niter=1;
end; % if

if ~isfield(S,'weightttest'), 
    S.weightttest=[];
end; 

if ~isfield(S,'testbands'),
    S.testbands=[];
end;

if ~isfield(S,'gridpos'),
    if ~isfield(S,'gridstep');
    S.gridstep = spm_input('Grid step (mm):', '+1', 'r', '5');
    end; 
end; % if


if ~isfield(S,'detrend'),
    S.detrend=[];
    end;
    
 if isempty(S.detrend),
     disp('detrending data by default');
     S.detrend=1;
 end; % 
 
 if ~isfield(S,'hanning'),
     S.hanning=[];
 end; 
    
 if isempty(S.hanning),
     disp('windowing data by default');
     S.hanning=1;
 end; % 

if ~isfield(S,'logflag'),
    S.logflag=[];
end; % if

if isempty(S.logflag),
    S.logflag=0;
end; % if
%  

if ~isfield(S,'regpc'),
    S.regpc=[];
end; % if

if isempty(S.regpc),
    S.regpc=0;
end; % if

if ~isfield(S,'suffix_str'),
    S.suffix_str=[];
end; % if


    
if ~isfield(S,'bootstrap'),
    S.bootstrap=[];
else
    Nboot=S.bootstrap;
    end;
    
if isempty(S.bootstrap),
    S.bootstrap=0;
    Nboot=1;
    end;
    
if ~isfield(S,'components'),
    S.components=[];
end; % if

if isempty(S.components),
    compind=1;
else
    compind=S.components;
end; % if


try
    D = S.D;
catch
    D = spm_select(1, 'mat', 'Select MEEG mat file');
    S.D = D;
end

if ischar(D)
    try
        D = spm_eeg_load(D);
    catch
        error(sprintf('Trouble reading file %s', D));
    end
end

[ok, D] = check(D, 'sensfid');

if ~ok
    if check(D, 'basic')
        errordlg(['The requested file is not ready for source reconstruction.'...
            'Use prep to specify sensors and fiducials.']);
    else
        errordlg('The meeg file is corrupt or incomplete');
    end
    return
end

modality = spm_eeg_modality_ui(D, 1, 1);

chanind=strmatch(modality, D.chantype);
chanind = setdiff(chanind,D.badchannels);
 channel_labels = D.chanlabels(chanind)';
 %chan_ind = setdiff(D.meegchannels('MEG'),D.badchannels)  


 if ~isfield(D,'inv')
     errordlg('Need to set up a forward model before you start');
 end;
 
if isfield(S, 'refchan') && ~isempty(S.refchan)
    refchan = S.refchan;
else
    refchan = [];
end

%% ============ Find or prepare head model

if ~isfield(D, 'val')
    D.val = 1;
end

   if ~isfield(S,'filenamestr'),
       S.filenamestr=[];
   end;%


for m = 1:numel(D.inv{D.val}.forward)
    if strncmp(modality, D.inv{D.val}.forward(m).modality, 3)
        vol  = D.inv{D.val}.forward(m).vol;
        if isa(vol, 'char')
            vol = ft_read_vol(vol);
        end
        datareg  = D.inv{D.val}.datareg(m);
    end
end


def_colormap=colormap;
jetmap=colormap('jet');
colormap(def_colormap);


if strcmp('EEG', modality)    
    sens = datareg.sensors;
    
else
    sens = D.sensors('MEG');    
    
end






Xdesign  =S.design.X;
c=S.design.contrast; %% c is contrast eg [ 0 1 -1] compare columns 2,3 of X


    
try S.design.X(:,1)-S.design.Xtrials-S.design.Xstartlatencies;
    catch
    error('Design windows missepcified');
end;

%X0  = X - X*c*pinv(c); 

%X0  = spm_svd(X0); %% X0 is null space i.e. everything that is happening in other columns of X


outfilenames='';


freqbands=[];
if ~isfield(S, 'freqbands')
    error('need to supply frequency bands')
end


if isempty(S.testbands),
    S.testbands=S.freqbands; %% bands to do the test on
end; % if
    
    
Nbands=numel(S.freqbands);




%% Check latencies are same here

%% now read in the first trial of data just to get sizes of variables right

%% now read in the first trial of data just to get sizes of variables right
Ntrials=size(S.design.X,1);
Isamples = D.indsample([S.design.Xstartlatencies(1) S.design.Xstartlatencies(1)+S.design.Xwindowduration]);
Nsamples= diff(Isamples)+1;
Nchans=length(channel_labels);

if S.hanning,
    fftwindow=hamming(Nsamples);
else
    disp('not windowing');
    fftwindow=ones(Nsamples,1);
end; 

allfftwindow=repmat(fftwindow,1,Nchans);
NumUniquePts = ceil((Nsamples+1)/2); %% data is real so fft is symmetric

if NumUniquePts<=2,
    error('Need to have more than 2 samples of data');
end;
fftnewdata=zeros(Ntrials,NumUniquePts,Nchans);
allepochdata=zeros(Ntrials,Nchans,Nsamples); %% for loading in data quickly

fHz = (0:NumUniquePts-1)*D.fsample/Nsamples;


if ~isfield(S,'Nfeatures'),
    Nfeatures=floor(Ntrials/3);
else
    Nfeatures=S.Nfeatures;
end;


%% now read in all trialtype and hold them as windowed fourier transforms
[uniquewindows]=unique(S.design.Xstartlatencies);
Nwindows=length(uniquewindows);
    

%% GET DATA- put each trial in allepochdata in same order as design matrix (i.e. remove dependence on Xtrials and Xstartlatencies)
TtofT=1; % e15; %% tesla to femto tesla - turned off
%%disp('rescaling from tesla to fT !!');
for i=1:Nwindows,     %% puts trials into epoch data according to order of design.X structures
    Isamples = D.indsample([uniquewindows(i) uniquewindows(i)+S.design.Xwindowduration]);
     useind=find(uniquewindows(i)==S.design.Xstartlatencies);
    Itrials =S.design.Xtrials(useind); %% indices into design.X structures
    allepochdata(useind,:,:)=permute(TtofT.*squeeze(D(D.indchannel(channel_labels), Isamples(1):Isamples(2), Itrials)), [3 1 2]); %% get an epoch of data with channels in columns
end; % for i

% for i=1:Nwindows,     %% puts trials into epoch data according to order of design.X structures
%     winstart=uniquewindows(i); %% window start
%     cfg=[];
%     cfg.keeptrials='yes';
%     cfg.channel=channel_labels;
%     cfg.feedback='off';
%     useind=find(winstart==S.design.Xstartlatencies); %% indices into design.X structures
%     cfg.trials=S.design.Xtrials(useind); %% trials starting at these times
%     cfg.latency=[winstart winstart+S.design.Xwindowduration];
%     subdata=ft_timelockanalysis(cfg,data); % subset of original data
%     if length(subdata.time)~=Nsamples,
%         error('Check the window specified is within epoch');
%         
%         end;
%     allepochdata(useind,:,:)=squeeze(subdata.trial); %% get an epoch of data with channels in columns
% end; % for i


for i=1:Ntrials, %% read in all individual trial types
    epochdata=squeeze(allepochdata(i,:,:))'; %% get an epoch of data with channels in columns
    if S.detrend==1
        dtepochdata=detrend(epochdata); %% detrend epoch data, this includes removind dc level. NB. This will have an effect on specific evoked response components !
        else 
        dtepochdata=epochdata; %% no dc removal, no detrend : this will have effect on accuracy of fourier estimate at non dc bins
        end; % detrend 
    wdtepochfft=dtepochdata.*allfftwindow; %% windowed 
    epochfft=fft(wdtepochfft);
    fftnewdata(i,:,:)=epochfft(1:NumUniquePts,:); % .*filtervect';    
end; 
if size(fftnewdata,3)~=Nchans,
    size(fftnewdata)
    error('Data dimension mismatch');
end;

if isfield(S,'return_data')
    stats.allepochdata=allepochdata;
    end;
    
clear allepochdata; %% no longer needed



%% now have an fft for each channel in each condition

    
% %%
cfg                       = [];

if ismember(modality, {'MEG', 'MEGPLANAR'})
    disp('Reducing possible source orientations to a tangential plane for MEG');
    cfg.reducerank = 2;
end;
cfg.grad=sens;
cfg.channel = channel_labels;
cfg.vol                   = vol;


 
if  isempty(S.gridpos),
    cfg.resolution            = S.gridstep;
else
    disp('USING pre-specified gridpoints');
    cfg.grid.pos=S.gridpos; %% predefined grid
    cfg.grid.inside=[1:size(S.gridpos,1)]; %% assume all in head
    cfg.grid.outside=[];
    end;

    
cfg.feedback='off';
cfg.inwardshift           = 0; % mm

if ~isfield(S,'grid'),
    disp('preparing leadfield');
    grid                      = ft_prepare_leadfield(cfg);
else
    disp('Using precomputed leadfield');
    grid=S.grid;
end; % if 


maskedgrid_inside_ind=[1:length(grid.inside)];
if ~isempty(S.maskgrid),
    if length(S.maskgrid)~=length(grid.inside),
        error('mask and grid points must be of same size');
        end;
    maskedgrid_inside_ind=find(S.maskgrid==1); %% indices into grid.inside 
    end; 

 if cfg.reducerank, %% follow up rank reduction and remove redundant dimension from lead fields
    for i=1:length(maskedgrid_inside_ind), %% 81
        lf1=cell2mat(grid.leadfield(grid.inside(maskedgrid_inside_ind(i))));
        [u1,s1,v1]=svd(lf1'*lf1);
        grid.leadfield(grid.inside(maskedgrid_inside_ind(i)))={lf1*u1(:,1:cfg.reducerank)};
        %normlf(i)=std(dot(lfnew',lfnew'));
    end;
 end; % if reduce rank
 
 %[a,b]=min(normlf) 
 %origin=grid.pos(grid.inside(maskedgrid_inside_ind(b)),:)

    
%% Now have all lead fields and all data
%% Now do actual beamforming
%% decide on the covariance matrix we need
%% construct covariance matrix within frequency range of interest

disp('now running through freq bands and constructing t stat images');
origfftnewdata=fftnewdata;
for boot=1:Nboot,
    
    bttrials=randi(Ntrials,Ntrials,1);
    if boot==1,
        bttrials=1:Ntrials;
    else
        disp(['Bootstrap run' num2str(boot)]);
    end; % if boot
    fftnewdata=origfftnewdata(bttrials,:,:);
    
for fband=1:Nbands,
    freqrange=S.freqbands{fband};
    freq_ind=intersect(find(fHz>=freqrange(1)),find(fHz<freqrange(2)));
    if length(freq_ind)<=1,
        disp(sprintf('Cannot use band %3.2f-%3.2f',freqrange(1),freqrange(2)));
        error('Need more than one frequency bin in the covariance band');
        end
    freqrangetest=S.testbands{fband};
    Ntestbands=length(freqrangetest)/2;
    if floor(Ntestbands)~=Ntestbands,
        error('Need pairs of values in testband');
        end;
    freq_indtest=[];
    freq_teststr='';
    for k=1:Ntestbands
        newtestind=intersect(find(fHz>=freqrangetest(k*2-1)),find(fHz<freqrangetest(k*2)));
        freq_indtest=[freq_indtest newtestind];    
        freq_teststr=sprintf('%s %3.2f-%3.2fHz,',freq_teststr,fHz(min(newtestind)),fHz(max(newtestind)));
        end; % for k
    if length(setdiff(freq_indtest,freq_ind))>0,
            error('Bands in test band are not within covariance band')
         end;
    covtrial=zeros(Nchans,Nchans);
   
    
    if ~isempty(S.weightttest),
        [fweightedt,weightspectindt]=intersect(fHz,S.weightttest(fband).fHz);
        if abs(max(fHz(freq_ind)-S.weightttest(fband).fHz))>0,
            error('weight ttest vector wrong length');
            end; % 
        tfiltervect=S.weightttest(fband).vect; %% weighted by previous mv analysis
        else
        tfiltervect=ones(length(freq_indtest),1);
        end;  % weightspect
        


for i=1:Ntrials, %% read in all individual trial types 
     ffttrial=squeeze(fftnewdata(i,freq_ind,:)); % .*Allfiltervect;
     covtrial=covtrial+real(cov(ffttrial));
end; % for i
  covtrial=covtrial/Ntrials;
  allsvd = svd(covtrial);
  cumpower=cumsum(allsvd)./sum(allsvd);
  nmodes99=min(find(cumpower>0.99));
  disp(sprintf('99 percent of the power in this data is in the first %d principal components of the cov matrix',nmodes99));
  disp(sprintf('largest/smallest eigenvalue=%3.2f',allsvd(1)/allsvd(end)));
  disp(sprintf('\nFrequency resolution %3.2fHz',mean(diff(fHz))));
  noise = allsvd(end); %% use smallest eigenvalue
  noise_id=eye(size(covtrial)).*noise;
  redNfeatures=Nfeatures;
  
    
  disp(sprintf('covariance band from %3.2f to %3.2fHz (%d bins), test band %s (%d bins)',fHz(freq_ind(1)),fHz(freq_ind(end)),length(freq_ind),freq_teststr,length(freq_indtest)))
  
    
  
  
  lambda = (S.regpc/100) * sum(allsvd)/size(covtrial,1); %% scale lambda relative mean eigenvalue
  disp(sprintf('regularisation =%3.2f percent',S.regpc));
  cinv=pinv(covtrial+eye(size(covtrial,1))*lambda); %% get inverse 
  
      



 tstat=zeros(length(grid.inside),S.Niter);
 normdiff=zeros(length(grid.inside),S.Niter);
maxt=zeros(2,S.Niter);
power_trial=zeros(Ntrials,length(freq_indtest));
evoked_trial=zeros(Ntrials,length(freq_indtest));

TrueIter=1; %% no permutation for this iteration
for j=1:S.Niter, %% set up permutations in advance- so perms across grid points are identical
    randind(j,:)=randperm(Ntrials);
    if j==TrueIter,
        randind(j,:)=1:Ntrials; % don't permute first run
        end;
    end;
      
 dfe=Ntrials-rank(Xdesign);  % df test
 

    
    for i=1:length(maskedgrid_inside_ind), %% 81
        lf=cell2mat(grid.leadfield(grid.inside(maskedgrid_inside_ind(i))));
        
        %% get optimal orientation- direct copy from Robert's beamformer_lcmv.m
        projpower_vect=pinv(lf'*cinv*lf);
        [u1,s1,v1]=svd(lf'*cinv*lf);
        usecomp=length(s1);
        
        [u, s, v] = svd(real(projpower_vect));
        eta = u(:,compind);
        
        lf  = lf * eta; %% now have got the lead field at this voxel, compute some contrast
        weights=lf'*cinv/(lf'*cinv*lf); %% CORRECT WEIGHTS CALC
        
        if S.return_weights
            stats(fband).ctf_weights(i,:)=weights;
        end
        
        for j=1:Ntrials, %% this has to be done at each location
            fdata=squeeze(fftnewdata(j,freq_indtest,:));
            
            fdatatrial=fdata*weights';
            evoked_trial(j,:)=fdatatrial;
            if S.logflag,
                power_trial(j,:)=log(fdatatrial.*conj(fdatatrial));
            else
                power_trial(j,:)=fdatatrial.*conj(fdatatrial); %%
            end; % i 

        end; % for j
        
        
        power_flag=1; %% only look at power for now
            if power_flag,
                Yfull=power_trial; %% univariate test later so just take the mean
            else
                Yfull=evoked_trial;
            end; % if power_flag
        
              
       %% Now permute the rows of X if necessary
        for iter=1:S.Niter,
        
          
            
             X=Xdesign(randind(iter,:),:); %% randind(1,:)=1, i.e. unpermuted
             if boot>1, %% have to also shuffle design matrix with data in bootstrap
                tmp=X;
                X=tmp(bttrials,:);
                end;
            
            % Contrast
            
             
             Yvals=mean(Yfull')';

             B  = pinv(X)*Yvals;

% t statistic and significance test
            RSS   = sum((Yvals - X*B).^2);
            MRSS  = RSS / dfe;
            SE    = sqrt(MRSS*(c*pinv(X'*X)*c'));
            
            tstat(maskedgrid_inside_ind(i),iter)=c*B./SE;
            normdiff(maskedgrid_inside_ind(i),iter)=c*B/(weights*noise_id*weights'); %% maybe a factor missing here
            
    
        end; % for Niter
        
         
        
        
     
        if i/100==floor(i/100)
            disp(sprintf('done t stats for %3.2f percent of freq band %d of %d, log=%d',100*i/length(maskedgrid_inside_ind),fband,Nbands,S.logflag));
        end; % if
    
    
  
end; % for grid points


      
stats(fband).tstat=tstat;
 maxt=max(tstat(:,TrueIter));
 mint=min(tstat(:,TrueIter));
stats(fband).fHz=fHz;

dispthresh_uv=max(stats(fband).tstat)/2;
if S.Niter>1,
    %% get corrected p values to t
        allglobalmax=squeeze(max(abs(stats(fband).tstat(:,1:end))));
        [sortglobalmax,sortglobalmaxind]=sort(allglobalmax','descend');
        
        stats(fband).corrpmax_tstat=find(TrueIter==sortglobalmaxind)/length(sortglobalmaxind);
        stats(fband).thresh05globalmax_tstat=sortglobalmax(round(length(sortglobalmaxind)*5/100),:);
        dispthresh_uv=stats(fband).thresh05globalmax_tstat; % display only significant effects
end; % if
  

    mnipositions = spm_eeg_inv_transform_points(D.inv{D.val}.datareg.toMNI, grid.pos(grid.inside(maskedgrid_inside_ind),:));
    gridpositions=grid.pos(grid.inside(maskedgrid_inside_ind),:);


    
    sMRI = fullfile(spm('dir'), 'canonical', 'single_subj_T1.nii');
    
    
    csource=grid; %% only plot and write out unpermuted iteration
    csource.pow_tstat(csource.inside) = tstat(:,TrueIter);
    csource.pow_tstat(csource.outside)=0;
    csource.pos = spm_eeg_inv_transform_points(D.inv{D.val}.datareg.toMNI, csource.pos);
    
    zeromean_images=0; %% leave this off for now.
    if zeromean_images==1,
        imgmean=mean(normdiff(:,TrueIter));
        imgstd=mean(normdiff(:,TrueIter));
        disp(sprintf('Removing mean value %3.2f from normdiff image (std=%3.2f)!',imgmean,imgstd));
        normdiff(:,TrueIter)=normdiff(:,TrueIter)-imgmean;
        end;
     csource.normdiff(csource.inside) =normdiff(:,TrueIter);
    csource.normdiff(csource.outside)=0;
    
    
     % CTF positions inside head
    ctf_inside = spm_eeg_inv_transform_points(D.inv{D.val}.datareg.toMNI, csource.pos(csource.inside,:));
    
    
    if isempty(S.gridpos), %% only write images if they use whole volume
         
    
    cfg1 = [];
    cfg1.sourceunits   = 'mm';
    cfg1.parameter = 'pow_tstat';
    cfg1.downsample = 1;
    sourceint_pow_tstat = ft_sourceinterpolate(cfg1, csource, sMRI);
    
     cfg1 = [];
     cfg1.sourceunits   = 'mm';
     cfg1.parameter = 'normdiff';
     cfg1.downsample = 1;
     sourceint_normdiff= ft_sourceinterpolate(cfg1, csource, sMRI);
%     
    
    
    
    
    
    %% else %% write out the data sets
    disp('writing images');
    
    cfg = [];
    cfg.sourceunits   = 'mm';
    cfg.parameter = 'pow';
    cfg.downsample = 1;
    
     featurestr='';
    if S.bootstrap,
        featurestr=sprintf('%s_bt%03d_',featurestr,boot);
        end; % if
    % write t stat
    dirname='tstatBf_images';
    if S.logflag,
        dirname=[dirname '_log'];
    end;

    res = mkdir(D.path, dirname);
    outvol = spm_vol(sMRI);
    outvol.dt(1) = spm_type('float32');
    
        outvol.fname= fullfile(D.path, dirname, ['spmT_' spm_str_manip(D.fname, 'r') '_' num2str(S.freqbands{fband}(1)) '-' num2str(S.freqbands{fband}(2)) 'Hz' S.filenamestr S.suffix_str featurestr '.nii']);
        
        stats(fband).outfile_pow_tstat=outvol.fname;
        outvol = spm_create_vol(outvol);
        spm_write_vol(outvol, sourceint_pow_tstat.pow_tstat);
         
        jetmap=colormap('jet');
        if (isfield(S, 'preview') && S.preview)
            spm_check_registration(sMRI)
            prop=0.4;
            colourmap=jetmap;
            spm_orthviews('Addtruecolourimage',1,outvol.fname,colourmap,prop,maxt,mint);
            disp('Press any key to continue');
            pause;
        end; % if preview
        outvol.fname= fullfile(D.path, dirname, ['spmNdiff_' spm_str_manip(D.fname, 'r') '_' num2str(S.freqbands{fband}(1)) '-' num2str(S.freqbands{fband}(2)) 'Hz' S.filenamestr S.suffix_str featurestr '.nii']);
        
         stats(fband).outfile_normdiff=outvol.fname;
         outvol = spm_create_vol(outvol);
         spm_write_vol(outvol, sourceint_normdiff.normdiff);

    
    end; % if ~S.gridpos
    
end; % for fband=1:Nbands

end; % for boot   

bootlist= fullfile(D.path, dirname, ['bootlist_'  spm_str_manip(D.fname, 'r') '_' num2str(S.freqbands{fband}(1)) '-' num2str(S.freqbands{fband}(2)) 'Hz' featurestr '.mat']);
save(bootlist,'bttrials');

     
end % function


