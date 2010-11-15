function [stats,talpositions,gridpositions,grid,fftnewdata,alllf,allepochdata]=spm_eeg_ft_beamformer_cva(S)
% Computes power-based beamformer image
% FORMAT [outfilenames,ctf_inside,ctf_weights]=spm_eeg_ft_beamformer_cva (S)
%
% S         MEEG object where coregistration has been performed.
%
% 
% Outputs (1) normalised power, (2) t-stat and (3) multivarariate
% (Hotellings)images. Uses Sekihara eigenval approach to choose optimal
% direction.
%
% outfilenames       Output filenames (for 1,2,3)
%
%                    The following fields are returned if you set 
%                    S.return_weights=1:
%
% ctf_inside         CTF locations inside head
% ctf_weights        Corresponding beamformer weight matrices
%                   fftnewdata is fft of data in trial order
% _______________________________________________________________________
% Copyright (C) 2009 Institute of Neurology, UCL

% Gareth Barnes
% $Id: spm_eeg_ft_beamformer_cva.m 4120 2010-11-15 11:47:36Z gareth $

[Finter,Fgraph] = spm('FnUIsetup','Multivariate LCMV beamformer for power', 0);
%%

%% ============ Load SPM EEG file and verify consistency
if nargin == 0
    S = [];
end

try
    D = S.D;
catch
    D = spm_select(1, '\.mat$', 'Select EEG mat file');
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
%try
    channel_labels = D.chanlabels(strmatch(modality, D.chantype))';
%catch
%    warning('ASSUMING ALL CHANNELS ARE MEG CHANNELS');
%    channel_labels=D.chanlabels;
%end;

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
%D.inv{1}.forward(1).voltype

 try
     vol = D.inv{D.val}.forward.vol;
     datareg = D.inv{D.val}.datareg;
 catch
     D = spm_eeg_inv_mesh_ui(D, D.val, [], 1);
     D = spm_eeg_inv_datareg_ui(D, D.val);
     datareg = D.inv{D.val}.datareg;
 end
 
for m = 1:numel(D.inv{D.val}.forward)
    if strncmp(modality, D.inv{D.val}.forward(m).modality, 3)
        vol  = D.inv{D.val}.forward(m).vol;
        if isa(vol, 'char')
            vol = fileio_read_vol(vol);
        end
        datareg  = D.inv{D.val}.datareg(m); %%  NEED TO KNOW HOW TO HANDLE DATAREG
        disp('USING LAST FORWARD MODEL !');
    end
end



% Return beamformer weights
if ~isfield(S,'return_weights')
    ctf_weights=[];
    S.return_weights=0;
    alllf=[];
end

if ~isfield(S,'Niter')
    S.Niter=[];
end; % if 

if isempty(S.Niter),
    S.Niter=1;
end; % if


if ~isfield(S,'design'),
   error('Design matrix required');
end; % if

X=S.design.X;
c=S.design.contrast; %% c is contrast eg [ 0 1 -1] compare columns 2,3 of X

if ~isfield(S,'ttest'),
    S.ttest=[];
    end; % if

if S.ttest,
    TWOSAMPLETEST=1
else
    TWOSAMPLETEST=0
    end;
    
try S.design.X(:,1)-S.design.Xtrials-S.design.Xstartlatencies;
    catch
    error('Design windows missepcified');
end;

X0  = X - X*c*pinv(c);  %% make sure X0 is orthogonal to X
Xdesign   = full(X*c);
X0  = spm_svd(X0); %% X0 is null space i.e. everything that is happening in other columns of X


outfilenames='';






freqbands=[];
if ~isfield(S, 'freqbands')
    for i = 1:spm_input('Number of frequency bands:', '+1', 'r', '1', 1)
        outstr=sprintf('Band %d [st end] in Hz ',i);
        S.freqbands{i} = spm_input(outstr, '+1', 'r', '', 2);
        freqbands =[freqbands;S.freqbands{i}'];
    end
else 
    freqbands=cell2mat(S.freqbands');
end

if ~isfield(S,'testbands'),
    S.testbands=[];
end;

if isempty(S.testbands),
    S.testbands=S.freqbands; %% bands to do the test on
end; % if
    
    
Nbands=numel(S.freqbands);

if ~isfield(S,'weightspect'), 
    S.weightspect=[];
end; 

if ~isfield(S,'weightttest'), 
    S.weightttest=[];
end; 




% if ~isempty(S.weightspect),
%         if Nbands>1,
%             error('cannot weight spectra for more than one freqband');
%         end; % if 
% end; % if

if ~isfield(S,'gridpos'),
    if ~isfield(S,'gridstep');
    S.gridstep = spm_input('Grid step (mm):', '+1', 'r', '5');
    end; 
end; % if

if ~isfield(S,'rankflag'),
    S.rankflag=[];
end; % if
if isempty(S.rankflag),
    S.rankflag=0;
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

%% READ IN JUST THE DATA WE need
%% ASSUME THAT INPUT IS AS FOLLOWS
%% a list of trial types and latency ranges for 2 conditions/periods (active and
%% rest for now)
%% each condition has an associated time window which must of equal sizes
%% SO will have
%% a list of trial indices and latencies, a list of trial types of the same
%% length

try
    data = D.ftraw; %% convert to field trip format- direct memory map 
catch
    disp('failed to read data directly.. going slowly');
   data = D.ftraw(0); %% convert to field trip format - file pointers
end; 
%% Check latencies are same here

%% now read in the first trial of data just to get sizes of variables right
cfg=[];
cfg.keeptrials='no';
cfg.trials=1; 


Ntrials=size(S.design.X,1);
cfg.latency=[S.design.Xstartlatencies(1) S.design.Xstartlatencies(1)+S.design.Xwindowduration];

subdata=ft_timelockanalysis(cfg,data); 
Nsamples=length(subdata.time);
%channel_labels = D.chanlabels(D.meegchannels(modality));
Nchans=length(channel_labels);
dfe=Ntrials-rank(X);

if S.hanning,
    fftwindow=hamming(Nsamples);
else
    disp('not windowing');
    fftwindow=ones(Nsamples,1);
end; 
allfftwindow=repmat(fftwindow,1,Nchans);
NumUniquePts = ceil((Nsamples+1)/2); %% data is real so fft is symmetric

fftnewdata=zeros(Ntrials,NumUniquePts,Nchans);
allepochdata=zeros(Ntrials,Nchans,Nsamples); %% for loading in data quickly

%fftnewdata=zeros(Ntrials,Nsamples,Nchans);

fHz = (0:NumUniquePts-1)*D.fsample/Nsamples;

if ~isfield(S,'Nfeatures'),
     S.Nfeatures=[];
end;


if isempty(S.Nfeatures),
    Nfeatures=floor(Ntrials/3);
else
    Nfeatures=S.Nfeatures;
end;


%% now read in all trialtype and hold them as windowed fourier transforms
[uniquewindows]=unique(S.design.Xstartlatencies);
Nwindows=length(uniquewindows);
    
    
for i=1:Nwindows,     %% puts trials into epoch data according to order of design.X structures
    winstart=uniquewindows(i); %% window start
    cfg=[];
    cfg.keeptrials='yes';
    cfg.channel=channel_labels;
    cfg.feedback='off';
    useind=find(winstart==S.design.Xstartlatencies); %% indices into design.X structures
    cfg.trials=S.design.Xtrials(useind); %% trials starting at these times
    cfg.latency=[winstart winstart+S.design.Xwindowduration];
    subdata=ft_timelockanalysis(cfg,data); % subset of original data
    allepochdata(useind,:,:)=squeeze(subdata.trial); %% get an epoch of data with channels in columns
end; % for i
% 
% for i=1:Nwindows,
%     winstart=uniquewindows(i); %% window start
%     cfg=[];
%     cfg.keeptrials='yes';
%     cfg.channel=channel_labels;
%     cfg.feedback='off';
%     cfg.trials=S.design.Xtrials(find(winstart==S.design.Xstartlatencies)); %% trials starting at these times
%     cfg.latency=[winstart winstart+S.design.Xwindowduration];
%     subdata=ft_timelockanalysis(cfg,data); % subset of original data
%     allepochdata(cfg.trials,:,:)=squeeze(subdata.trial); %% get an epoch of data with channels in columns
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
% 
% for i=1:Ntrials, %% read in all individual trial types
% %     cfg.trials=S.design.Xtrials(i);
% %     cfg.latency=[S.design.Xstartlatencies(i) S.design.Xstartlatencies(i)+S.design.Xwindowduration]
% % 
% % 
% %     cfg.channel=channel_labels;
% %     cfg.feedback='off';
% %     
% %     %%    cfg.trl=cfg.trials;
% %     subdata=ft_timelockanalysis(cfg,data); % subset of original data
%     epochdata=squeeze(allepochdata(S.design.Xtrials(i),:,:))'; %% get an epoch of data with channels in columns
%     if S.detrend==1
%         dtepochdata=detrend(epochdata); %% detrend epoch data, this includes removind dc level. NB. This will have an effect on specific evoked response components !
%         else 
%         dtepochdata=epochdata; %% no dc removal, no detrend : this will have effect on accuracy of fourier estimate at non dc bins
%         end; % detrend 
%     wdtepochfft=dtepochdata.*allfftwindow; %% windowed
%     
%     epochfft=fft(wdtepochfft);
%     %%     epochcov=cov(epochfft);
%     
%     fftnewdata(i,:,:)=epochfft(1:NumUniquePts,:); % .*filtervect';
%    
%     
% end; 
% if size(fftnewdata,3)~=Nchans,
%     fftnewdata
%     error('Data dimension mismatch');
% end;

%clear allepochdata; %% no longer needed



%% now have an fft for each channel in each condition

    
% %%
cfg                       = [];
if strcmp('EEG', modality)
    cfg.elec = D.inv{D.val}.datareg.sensors;
    cfg.reducerank=3;
else
    cfg.grad = D.sensors('MEG');
    cfg.reducerank=2;
    disp('Reducing possible source orientations to a tangential plane for MEG');
end


cfg.channel = channel_labels;
cfg.vol                   = vol;

if ~isfield(S,'gridpos'),
    S.gridpos=[];
    end;
if ~isfield(S,'gridori'),
    S.gridori=[];
end;

 if ~isfield(S,'maskgrid'),
    S.maskgrid=[];
    end;

    
 if ~isfield(S,'dimauto'),
    S.dimauto=[];
    end;
 
    if ~isfield(S,'maskfile'),
        S.maskfile=[];
    end;
    
    if ~isfield(S,'randomweight'),
        S.randomweight=0;
    end;

 if ~isfield(S,'fixedweights'),
        S.fixedweights=[];
    end;

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

 if ~isempty(S.maskfile),
     if ~isempty(S.maskgrid),
         error('cannot have two masks defined');
        end;
     alltalpositions = spm_eeg_inv_transform_points(datareg.toMNI, grid.pos(grid.inside,:));
      V_aal = spm_vol(S.maskfile);
     [Y_aal,XYZ]=spm_read_vols(V_aal);
      mask_coords=XYZ(:,find(Y_aal>0));
      
    
    %% now express mni grid and mesh grid at same coarse scale to get
    %% intersection
      coarse_mm=10;
      mask_coarse=unique(round(mask_coords'/coarse_mm)*coarse_mm,'rows');
      grid_coarse= round(alltalpositions/coarse_mm)*coarse_mm;
      [overlap_pos,maskedgrid_inside_ind]=intersect(grid_coarse,mask_coarse,'rows'); %% mesh_vert_ind are mesh points in this mask
      maskedgrid_inside_ind=sort(maskedgrid_inside_ind);
 end; % if 
%  if isfield(S,'gridstep_Ylim'),
%     maskedgrid_inside_ind=find(positio
%     end;

if ~isempty(S.fixedweights),
    disp('NB USING SUPPLIED FIXED WEIGHTS');
    if size(S.fixedweights,1)~=length(maskedgrid_inside_ind),
        error('Suppied weights are the wrong size for grid');
        end;
end;


%% Now have all lead fields and all data
%% Now do actual beamforming
%% decide on the covariance matrix we need
%% construct covariance matrix within frequency range of interest

disp('now running through freq bands and constructing t stat images');
for fband=1:Nbands,
    freqrange=freqbands(fband,:);
    freq_ind=intersect(find(fHz>=freqrange(1)),find(fHz<freqrange(2)));
    freqrangetest=cell2mat(S.testbands(fband));
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
   
    if ~isempty(S.weightspect),
        [fweighted,weightspectind]=intersect(fHz,S.weightspect(fband).fHz);
        if abs(max(fHz(freq_ind)-S.weightspect(fband).fHz))>0,
            error('weight spect vector wrong length');
            end; % 
        filtervect=S.weightspect(fband).fAmp; %% arb filter
    else
        filtervect=ones(length(freq_ind),1);
    end;  % weightspect
    
    if ~isempty(S.weightttest),
        [fweightedt,weightspectindt]=intersect(fHz,S.weightttest(fband).fHz);
        if abs(max(fHz(freq_ind)-S.weightttest(fband).fHz))>0,
            error('weight ttest vector wrong length');
            end; % 
        tfiltervect=S.weightttest(fband).vect; %% weighted by previous mv analysis
        else
        tfiltervect=ones(length(freq_indtest),1);
        end;  % weightspect
        
    Allfiltervect=repmat(filtervect,1,Nchans);

for i=1:Ntrials, %% read in all individual trial types 
     ffttrial=squeeze(fftnewdata(i,freq_ind,:)).*Allfiltervect;
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
  redNfeatures=Nfeatures;
  
  if S.dimauto, %% only do this flag set
    if nmodes99<redNfeatures,
        redNfeatures=nmodes99+1;
        disp(sprintf('reducing number of features to those describing 99 percent of data (from %d to %d)',Nfeatures,redNfeatures));
    end; % 
   end; % if S.dimauto
    
    if length(freq_indtest)<redNfeatures,    
        disp(sprintf('reducing number of features to match bandwidth (from %d to %d) !!',redNfeatures,length(freq_indtest)));
        redNfeatures=length(freq_indtest);
    end;
  
    
  disp(sprintf('covariance band from %3.2f to %3.2fHz (%d bins), test band %s (%d bins)',fHz(freq_ind(1)),fHz(freq_ind(end)),length(freq_ind),freq_teststr,length(freq_indtest)))
  disp(sprintf('Using %d features',redNfeatures));
    
  
  
  lambda = (S.regpc/100) * sum(allsvd)/size(covtrial,1); %% scale lambda relative mean eigenvalue
  disp(sprintf('regularisation =%3.2f percent',S.regpc));
  cinv=pinv(covtrial+eye(size(covtrial,1))*lambda); %% get inverse - if N features has been reduced should maybe think of sorting this out too.
  
      

CVA_maxstat=zeros(length(grid.inside),2,S.Niter);
roymax=zeros(size(CVA_maxstat));

 tstat=zeros(length(grid.inside),2,S.Niter);
maxcva=zeros(2,S.Niter);
power_trial=zeros(Ntrials,length(freq_indtest));
evoked_trial=zeros(Ntrials,length(freq_indtest));

TrueIter=1; %% no permutation for this iteration
for j=1:S.Niter, %% set up permutations in advance- so perms across grid points are identical
    randind(j,:)=randperm(Ntrials);
    if j==TrueIter,
        randind(j,:)=1:Ntrials; % don't permute first run
        end;
    end;
    
  
  
  
  
    for i=1:length(maskedgrid_inside_ind), %% 81
        %disp('gridpoint')
        %i
        %length(maskedgrid_inside_ind)
        lf=cell2mat(grid.leadfield(grid.inside(maskedgrid_inside_ind(i))));
        
        %% get optimal orientation- direct copy from Robert's beamformer_lcmv.m
        projpower_vect=pinv(lf'*cinv*lf);
        if isempty(S.gridori),
            [u, s, v] = svd(real(projpower_vect)); %% MAY NOT BE POWER WE NEED TO OPTIMISE
            eta = u(:,1);
        else
            eta=S.gridori(grid.inside(maskedgrid_inside_ind(i)),:)';
        end;
%         if i==80,
%             disp('FIXING ORIENTATION AT SOURCE LOCATION !!')
%             eta=[0.9513    0.1621         0]';
%             end; % if i
        lf  = lf * eta; %% now have got the lead field at this voxel, compute some contrast
       
        weights=(lf'*cinv*lf)\lf'*cinv; %% no regularisation for now
        
        if S.fixedweights,
            weights=S.fixedweights(i,:);
            end;
        if S.randomweight,
            disp('RANDOM WEIGHT !');
            weights=randn(size(weights));
        end;
        %if S.return_weights
            stats(fband).ctf_weights(maskedgrid_inside_ind(i),:)=weights;
            alllf(maskedgrid_inside_ind(i),:)=lf;
        %end
        
        for j=1:Ntrials, %% this non-linear step (power estimation) has to be done at each location
            fdata=squeeze(fftnewdata(j,freq_indtest,:));
              if length(freq_indtest)==1,
                  fdata=fdata';
                  end; % if length
             fdatatrial=fdata*weights';
            evoked_trial(j,:)=fdatatrial;
            if S.logflag,
                power_trial(j,:)=log(fdatatrial.*conj(fdatatrial));
            else
                power_trial(j,:)=fdatatrial.*conj(fdatatrial); %%
            end; % i 

        end; % for j
        
        

        for power_flag=0:1,
            if power_flag,
                Yfull=power_trial;
            else
            %    Yfull=power_trial./repmat(std(power_trial),size(power_trial,1),1); 
                Yfull=[real(evoked_trial), imag(evoked_trial)]; %% split into sin and cos parts 
            end; % if power_flag
            %Yfull=power_trial;
        Y     = Yfull - X0*(X0'*Yfull); %% eg remove DC level or drift terms from all of Y
        [u s v] = spm_svd(Y,-1,-1); % ,1/4);              %% get largest variance in Y   -1s make sure nothing is removed
            %% reduce dimensions to Nfeatures
        thresh=s(redNfeatures,redNfeatures);
        if redNfeatures==length(freq_indtest),
            thresh=0;
            end;
      
          [u,s,v]=svd(Y'*Y); %% reduce dimensions of useful signal
       %   U=u(1:redNfeatur es,:);
          U=u(:,1:redNfeatures); %% IMPORTANT ?
          Y   = Y*U;
        
  %% NOTE TSTAT COMPUTED Y BEFORE DIMENSION REDUCTION- COMPARISON ONLY VALID WHEN Y
  %% and YFULL have same dimensionality
  %% compare with hotellings
        
       %% Now permute the rows of X if necessary
        for iter=1:S.Niter,
        
            
        X=Xdesign(randind(iter,:),:); %% randind(1,:)=1, i.e. unpermuted
  
              %-Canonical Variates Analysis
    %   ==========================================================================
    % remove null space of contrast
     %--------------------------------------------------------------------------
            Y     = Y - X0*(X0'*Y); %% eg remove DC level or drift terms from all of Y
            X     = X - X0*(X0'*X);
            P     = pinv(X);
  
        if TWOSAMPLETEST,
            Ym=mean(Y,2);
            B  = pinv(X)*Ym;
            RSS   = sum((Ym - X*B).^2);
            MRSS  = RSS / dfe;
            SE    = sqrt(MRSS*(pinv(X'*X)));
            tstat(maskedgrid_inside_ind(i),power_flag+1,iter)=B./SE;
            Ftstat=(B./SE).^2;
        end; % if TWOSAMPLETEST
    
        
        [n,b] = size(X);                       
        [n,m] = size(Y); %% n is number of epochs, m is number of features
        b     = rank(X);
        h     = min(b,m); %% either number of features or rank of X
        f     = n - b - size(X0,2); %% number of epochs- rank(X) - number of confounds
 
 
        % generalised eigensolution for treatment and residual sum of squares
        %--------------------------------------------------------------------------
        T     = X*(P*Y); %% projection of X onto Y ?
        SST   = T'*T;
        SSR   = Y - T; %% residuals in Y (unexplained by X)
        SSR   = SSR'*SSR;
        [v,d] = eig(SSR\SST); %% explained/unexplained variance
        [q,r] = sort(-real(diag(d)));
        r     = r(1:h);
        d     = real(d(r,r));
        v     = real(v(:,r));
        V     = U*v;                          % canonical vectors  (data)
        v     = Y*v;                          % canonical variates (data)
        W     = P*v;                          % canonical vectors  (design)
        w     = X*W;                          % canonical variates (design)
        C     = c*W;                          % canonical contrast (design)
        
        
        
 
% inference on dimensionality - p(i) test of D >= i; Wilk's Lambda := p(1)
%--------------------------------------------------------------------------
        roymax(maskedgrid_inside_ind(i),power_flag+1,iter)=d(1,1); %% maximum root is the first one
        Roy2F=(f+m-max(b,m))/max(b,m); %% conversion from Roy's max root to F
        

        cval  = log(diag(d) + 1);
        for i1 = 1:h
            chi(i1) = (f - (m - b + 1)/2)*sum(cval(i1:h));
            df(i1)  = (m - i1 + 1)*(b - i1 + 1); % m is the number of features, b is the rank of the design matrix
            p(i1)   = 1 - spm_Xcdf(chi(i1),df(i1));
            p05thresh(i1) = spm_invXcdf(1-0.05,df(i1));
        end
        
       
        CVA_maxstat(maskedgrid_inside_ind(i),power_flag+1,iter)=chi(1); %% get wilk's lambda stat
        CVA_otherdim(maskedgrid_inside_ind(i),power_flag+1,iter,1:h-1)=chi(2:h); %% tests on other dimensions
        
        allw(maskedgrid_inside_ind(i),power_flag+1,iter,:,:)=W;
        
        
        if CVA_maxstat(maskedgrid_inside_ind(i),power_flag+1,iter)>maxcva(power_flag+1,iter),
            stats(fband).bestchi(power_flag+1,1:h,iter)=chi;
            if length(V)>1,
                if power_flag,
                 stats(fband).bestVpw(1:h,iter,:)=V';
                 stats(fband).bestvpw(1:h,iter,:)=v';
                 stats(fband).bestUpw=U;
                else
                 stats(fband).bestVev(1:h,iter,:)=V';
                 stats(fband).bestvev(1:h,iter,:)=v';
                 stats(fband).bestUev=U;
                end;
            end;% 
           
            stats(fband).pval(power_flag+1,1:h,iter)=p;
            stats(fband).df(power_flag+1,1:h,iter)=df;
            stats(fband).p05alyt(power_flag+1,1:h)=p05thresh;
            stats(fband).maxind(power_flag+1,iter)=maskedgrid_inside_ind(i);
            stats(fband).freqHz=fHz(freq_indtest);
            stats(fband).freq_indtest=freq_indtest;
            stats(fband).freq_ind=freq_ind;
            stats(fband).maxw=w;
            stats(fband).maxW=W;
            maxcva(power_flag+1,iter)=CVA_maxstat(maskedgrid_inside_ind(i),power_flag+1,iter);
            end; % if CVA_Stat
        if (iter==TrueIter) && (length(V)>1),
            if power_flag,
                stats(fband).allVpw(1:h,i,:)=V';
                else
                stats(fband).allVev(1:h,i,:)=V';
                end;
            end;% if iter
        
        end; % for Niter
        
         
        
        end; % for power_flag  
     
        if i/100==floor(i/100)
            disp(sprintf('done CVA stats for %3.2f percent of freq band %d of %d, log=%d, rank=%d',100*i/length(maskedgrid_inside_ind),fband,Nbands,S.logflag,S.rankflag));
        end; % if
    
    
  
end; % for grid points


      
stats(fband).maskedgrid_inside_ind=maskedgrid_inside_ind;
stats(fband).CVAmax=CVA_maxstat;
stats(fband).CVAotherdim=CVA_otherdim;
stats(fband).roymax=roymax;
stats(fband).Roy2F=Roy2F;
stats(fband).allw=allw;
stats(fband).tstat=tstat.^2;
stats(fband).fHz=fHz;
stats(fband).dfmisc=[b h m];
dispthresh_mv=max(stats(fband).CVAmax)/2; % default display thresholds
max_mv=max(stats(fband).CVAmax); % default display thresholds
dispthresh_uv=max(stats(fband).tstat)/2;
max_uv=max(stats(fband).CVAmax); % default display thresholds
if S.Niter>1,
    %% get corrected p values for CVA
    allglobalmax=squeeze(max(stats(fband).CVAmax(:,:,1:end)));
    [sortglobalmax,sortglobalmaxind]=sort(allglobalmax','descend');
    allglobalmax_roy=squeeze(max(stats(fband).roymax(:,:,1:end)));
    [sortglobalmax_roy,sortglobalmaxind_roy]=sort(allglobalmax_roy','descend');
    
    for k=1:2,
        stats(fband).corrpmax_cva(k)=find(sortglobalmaxind(:,k)==TrueIter)/length(sortglobalmaxind);
        stats(fband).corrpmax_roy(k)=find(sortglobalmaxind_roy(:,k)==TrueIter)/length(sortglobalmaxind_roy);
        end;
     stats(fband).thresh05globalmax_cva=sortglobalmax(round(length(sortglobalmaxind)*5/100),:);
     stats(fband).thresh05globalmax_roy=sortglobalmax_roy(round(length(sortglobalmaxind)*5/100),:);
    dispthresh_mv=stats(fband).thresh05globalmax_cva; % display only significant effects
    if TWOSAMPLETEST
    %% get corrected p values to t
        allglobalmax=squeeze(max(abs(stats(fband).tstat(:,:,1:end))));
        [sortglobalmax,sortglobalmaxind]=sort(allglobalmax','descend');
        for k=1:2,
            stats(fband).corrpmax_tstat(k)=find(sortglobalmaxind(:,k)==TrueIter)/length(sortglobalmaxind); 
            end;
        stats(fband).thresh05globalmax_tstat=sortglobalmax(round(length(sortglobalmaxind)*5/100),:);
        dispthresh_uv=stats(fband).thresh05globalmax_tstat; % display only significant effects
    end; % if TWOSAMPLETEST
end; % if
  
%save('tmp_data_cva.mat');
talpositions = spm_eeg_inv_transform_points(datareg.toMNI, grid.pos(grid.inside,:)); %% MADE DATAREG STAND ALONE
gridpositions=grid.pos(grid.inside,:);


    
    sMRI = fullfile(spm('dir'), 'canonical', 'single_subj_T1.nii');
    
    
    csource=grid; %% only plot and write out unpermuted iteration
    csource.pow_maxchi(csource.inside) = CVA_maxstat(:,2,TrueIter);  %power
    csource.evoked_maxchi(csource.inside) = CVA_maxstat(:,1,TrueIter);
    csource.pow_tstat(csource.inside) = tstat(:,2,TrueIter).^2;
    
    csource.pow_maxchi(csource.outside)=0;
    csource.evoked_maxchi(csource.outside)=0;
    csource.pow_tstat(csource.outside)=0;
    
    csource.pos = spm_eeg_inv_transform_points(datareg.toMNI, csource.pos);
    
    
    
    if isempty(S.gridpos), %% only write images if they use whole volume
        
    
    cfg1 = [];
    cfg1.sourceunits   = 'mm';
    cfg1.parameter = 'pow_maxchi';
    cfg1.downsample = 1;
    sourceint_pow_maxchi = ft_sourceinterpolate(cfg1, csource, sMRI);
    
    cfg1 = [];
    cfg1.sourceunits   = 'mm';
    cfg1.parameter = 'evoked_maxchi';
    cfg1.downsample = 1;
    sourceint_evoked_maxchi = ft_sourceinterpolate(cfg1, csource, sMRI);
    
    if TWOSAMPLETEST
    cfg1 = [];
    cfg1.sourceunits   = 'mm';
    cfg1.parameter = 'pow_tstat';
    cfg1.downsample = 1;
    sourceint_pow_tstat = ft_sourceinterpolate(cfg1, csource, sMRI);
    end; % IF TWOSAMPLETEST
    
    
    %[maxvals,maxind]=max(stats(fband).ctf_weights');
    % [minvals,minind]=max(-stats(fband).ctf_weights');
    [maxvals,maxind]=max(alllf');
    [minvals,minind]=max(-alllf');
     
       weight_extrema=[maxind' minind'];
       weight_extrema=sort(weight_extrema')'; %% doesn't matter if extrema are reversed
       u_extrema=unique(weight_extrema,'rows');
       stats(fband).Nu_extrema=size(u_extrema,1);
       alpha=0.05;
       stats(fband).alyt_thresh_Chi_05=spm_invXcdf(1-alpha/(stats(fband).Nu_extrema),Nfeatures*b); %% estimate of volumetric threshold
    
    if isfield(S,'S.dispthresh_chi'),
        dispthresh_mv(2)=S.dispthresh_chi;
        end;
    
        
    %end; % if preview
    
    
    %% else %% write out the data sets
    disp('writing images');
    
    cfg = [];
    cfg.sourceunits   = 'mm';
    cfg.parameter = 'pow';
    cfg.downsample = 1;
    % write t stat
    dirname='cvaBf_images';
    if S.logflag,
        dirname=[dirname '_log'];
    end;

    res = mkdir(D.path, dirname);
    outvol = spm_vol(sMRI);
    outvol.dt(1) = spm_type('float32');
    featurestr=[S.filenamestr 'Nf' num2str(redNfeatures)] ;
    outvol.fname= fullfile(D.path, dirname, ['chi_pw_'  spm_str_manip(D.fname, 'r') '_' num2str(freqbands(fband,1)) '-' num2str(freqbands(fband,2)) 'Hz' featurestr '.nii']);
    stats(fband).outfile_chi_pw=outvol.fname;
    outvol = spm_create_vol(outvol);
    spm_write_vol(outvol, sourceint_pow_maxchi.pow_maxchi);
    cmap=colormap;
    jetmap=colormap('jet');
    if (isfield(S, 'preview') && S.preview)
            spm_check_registration(sMRI)
            prop=0.4;
            colourmap=jetmap;
            spm_orthviews('Addtruecolourimage',1,outvol.fname,colourmap,prop,max_mv(2),dispthresh_mv(2));
            disp('Chi pw image. Press any key to continue..');
            pause;
        end; % if preview
    colormap(cmap);
    outvol.fname= fullfile(D.path, dirname, ['chi_ev_' spm_str_manip(D.fname, 'r') '_' num2str(freqbands(fband,1)) '-' num2str(freqbands(fband,2)) 'Hz' featurestr '.nii']);
    stats(fband).outfile_chi_ev=outvol.fname;
    outvol = spm_create_vol(outvol);
    spm_write_vol(outvol, sourceint_evoked_maxchi.evoked_maxchi);
    
        if (isfield(S, 'preview') && S.preview)
            spm_check_registration(sMRI)
            prop=0.4;
            colourmap=jetmap;
            spm_orthviews('Addtruecolourimage',1,outvol.fname,colourmap,prop,max_mv(1),dispthresh_mv(1));
            disp('Chi ev image. Press any key to continue..');
            pause;
        end; % if preview
        
    colormap(cmap);
    if TWOSAMPLETEST
        outvol.fname= fullfile(D.path, dirname, ['tstat_pow_' spm_str_manip(D.fname, 'r') '_' num2str(freqbands(fband,1)) '-' num2str(freqbands(fband,2)) 'Hz' S.filenamestr '.nii']);
        stats(fband).outfile_pow_tstat=outvol.fname;
        outvol = spm_create_vol(outvol);
        spm_write_vol(outvol, sourceint_pow_tstat.pow_tstat);
    end; % if twosample
    
%     
    end; % if ~S.gridpos
    
    
       
    
end; % for fband=1:Nbands
    
end % function


