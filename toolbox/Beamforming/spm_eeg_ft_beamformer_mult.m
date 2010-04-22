function [outfilenames,ctf_inside,ctf_weights,fftnewdata]=spm_eeg_ft_beamformer_mult(S)
% Compute power-based beamformer image
% FORMAT [outfilenames,ctf_inside,ctf_weights]=spm_eeg_ft_beamformer_mult(S)
%
% S            - struct (optional)
% (optional) fields of S:
% S.D          - meeg object or filename
%
% Outputs (1) normalised power, (2) t-stat and (3) multivarariate
% (Hotellings)images. Uses Sekihara eigenval approach to choose optimal
% direction.
%
% outfilenames - Output filenames (for 1,2,3)
%
%                The following fields are returned if you set 
%                S.return_weights=1:
%
% ctf_inside   - CTF locations inside head
% ctf_weights  - Corresponding beamformer weight matrices
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: spm_eeg_ft_beamformer_mult.m 3833 2010-04-22 14:49:48Z vladimir $

[Finter,Fgraph] = spm('FnUIsetup','Multivariate LCMV beamformer for power', 0);
%%

%% ============ Load SPM EEG file and verify consistency
if nargin == 0
    S = [];
end

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
channel = D.chanlabels(strmatch(modality, D.chantype))';

if isfield(S, 'refchan') && ~isempty(S.refchan)
    refchan = S.refchan;
else
    refchan = [];
end

%% ============ Find or prepare head model

if ~isfield(D, 'val')
    D.val = 1;
end

D.inv{1}.forward(1).voltype

for m = 1:numel(D.inv{D.val}.forward)
    if strncmp(modality, D.inv{D.val}.forward(m).modality, 3)
        vol  = D.inv{D.val}.forward(m).vol;
        if isa(vol, 'char')
            vol = ft_read_vol(vol);
        end
        datareg  = D.inv{D.val}.datareg(m);
    end
end

 try
     vol = D.inv{D.val}.forward.vol;
     datareg = D.inv{D.val}.datareg;
 catch
     D = spm_eeg_inv_mesh_ui(D, D.val, [], 1);
     D = spm_eeg_inv_datareg_ui(D, D.val);
     datareg = D.inv{D.val}.datareg;
 end

% Return beamformer weights
if ~isfield(S,'return_weights')
    ctf_weights=[];
    S.return_weights=0;
end
 

clb = D.condlist;
triallist=[];
trialtypes=[];
latencies=[];
outfilenames='';

contrast_str='';
if ~isfield(S,'timewindows'),
    S.timewindows={[],[]};
end; % if

Ntimewindows=2; % %% fixed for now
for i = 1:Ntimewindows, %%spm_input('Number of time windows:', '+1', 'r', '1', 1)
    
    contrast_str=[contrast_str 'Tr']
    if isempty(S.timewindows{i});
        [selection, ok]= listdlg('ListString', clb, 'SelectionMode', 'multiple' ,'Name', 'Select conditions' , 'ListSize', [400 300]);
        if ~ok
            return;
        end
        S.trigger{i}=clb(selection);
        outstr=sprintf('Offset (sec) from %s',cell2mat(S.trigger{i}));
        starttime= spm_input(outstr, '+1', 'r', '', 1);
        if i==1,
            outstr=sprintf('Duration (sec) ');
            duration= spm_input(outstr, '+1', 'r', '', 1);
        end; % if i==1
        S.timewindows{i}=[starttime;starttime+duration];
        S.preview=0;
    end; % if isfield timewindows
    
    Ntrialspercond=0;
    
    for s=1:numel(S.trigger{i}), %% maybe multiple triggers in selection
        trigname=cell2mat(S.trigger{i}(s));
        Ntrialspercond=Ntrialspercond+length(D.pickconditions(trigname)); %% number of trials for this condition
        contrast_str=[contrast_str trigname];
        triallist=[triallist , D.pickconditions(trigname)];
        contrast_str=[contrast_str ','];
    end; % for selection
    trialtypes=[trialtypes ; ones(Ntrialspercond,1).*i];
    latencies =[latencies;repmat([S.timewindows{i}],1,Ntrialspercond)'];
    contrast_str=[contrast_str num2str(S.timewindows{i}(1)) '-' num2str(S.timewindows{i}(2)) 'sec'];
    if i==1,
        contrast_str=[contrast_str '_vs'];
    end; % if i
end; % for i

contrast_str=sprintf('%s',deblank(contrast_str));

Nconditions=numel(S.timewindows);


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
Nbands=numel(S.freqbands);

if ~isfield(S,'gridstep');
 S.gridstep = spm_input('Grid step (mm):', '+1', 'r', '5');
end; 

if ~isfield(S,'rankflag'),
    S.rankflag=[];
end; % if
if isempty(S.rankflag),
    S.rankflag=0;
end; % if

if ~isfield(S,'logflag'),
    S.logflag=[];
end; % if
if isempty(S.logflag),
    S.logflag=0;
end; % if
%  
% if ~isfield(S, 'lambda')
%     S.lambda = [num2str(spm_input('Regularization (%):', '+1', 'r', '0')) '%'];
% end


%% READ IN JUST THE DATA WE need
%% ASSUME THAT INPUT IS AS FOLLOWS
%% a list of trial types and latency ranges for 2 conditions/periods (active and
%% rest for now)
%% each condition has an associated time window which must of equal sizes
%% SO will have
%% a list of trial indices and latencies, a list of trial types of the same
%% length

try
    data = D.ftraw(0); %% convert to field trip format- direct memory map 
catch
    disp('failed to read data directly.. going slowly');
   data = D.ftraw(1); %% convert to field trip format - file pointers
end; 
%% Check latencies are same here

%% now read in the first trial of data just to get sizes of variables right
cfg=[];
cfg.keeptrials='no';
cfg.trials=1; 
cfg.latency=latencies(1,:)
subdata=ft_timelockanalysis(cfg,data); 
Nsamples=length(subdata.time);
channel_labels = D.chanlabels(D.meegchannels(modality));
Nchans=length(channel_labels);
Ntrials=length(trialtypes);
fftwindow=hamming(Nsamples);
allfftwindow=repmat(fftwindow,1,Nchans);
NumUniquePts = ceil((Nsamples+1)/2); %% data is real so fft is symmetric
fftnewdata=zeros(Ntrials,NumUniquePts,Nchans);

fHz = (0:NumUniquePts-1)*D.fsample/Nsamples;



%% now read in all trialtype and hold them as windowed fourier transforms
cfg=[];
cfg.keeptrials='yes';
for i=1:length(trialtypes), %% read in all individual trial types
    cfg.trials=triallist(i);
    cfg.latency=latencies(i,:);
    cfg.channel=channel_labels;
    cfg.feedback='off';
    %%    cfg.trl=cfg.trials;
    subdata=ft_timelockanalysis(cfg,data); % subset of original data
    epochdata=squeeze(subdata.trial)'; %% get an epoch of data with channels in columns
    dtepochdata=detrend(epochdata); %% detrend epoch data, this includes removind dc level
    wdtepochfft=dtepochdata.*allfftwindow; %% windowed
    origcov=cov(wdtepochfft);
    epochfft=fft(wdtepochfft);
    %%     epochcov=cov(epochfft);
    fftnewdata(i,:,:)=epochfft(1:NumUniquePts,:);
    
end; 




%% now have an fft for each channel in each condition
    
% %%
cfg                       = [];
if strcmp('EEG', modality)
    cfg.elec = D.inv{D.val}.datareg.sensors;
    cfg.reducerank=3;
else
    cfg.grad = D.sensors('MEG');
    cfg.reducerank=2;
    disp('reducing possible source orientations to a tangential plane for MEG');
end

cfg.channel = channel_labels;
cfg.vol                   = vol;
cfg.resolution            = S.gridstep;
cfg.feedback='off';
cfg.inwardshift           = 0; % mm
grid                      = ft_prepare_leadfield(cfg);



%% Now have all lead fields and all data
%% Now do actual beamforming
%% decide on the covariance matrix we need
%% construct covariance matrix within frequency range of interest

disp('now running through freq bands and constructing t stat images');
for f=1:Nbands,
    freqrange=freqbands(f,:);
    freq_ind=intersect(find(fHz>=freqrange(1)),find(fHz<freqrange(2)));
    covtrial=zeros(Nchans,Nchans);
    
    for i=1:length(trialtypes), %% read in all individual trial types
        ffttrial=squeeze(fftnewdata(i,:,:));
        covtrial=covtrial+real(cov(ffttrial(freq_ind,:)));
    end; % for i
    covtrial=covtrial./length(trialtypes);
    allsvd = svd(covtrial);
    cumpower=cumsum(allsvd)./sum(allsvd);
    nmodes99=min(find(cumpower>0.99));
    disp(sprintf('99 percent of the power in this data is in the first %d principal components of the cov matrix',nmodes99));
    disp(sprintf('largest/smallest eigenvalue=%3.2f',allsvd(1)/allsvd(end)));
    disp(sprintf('\nFrequency resolution %3.2fHz',mean(diff(fHz))));
    noise = allsvd(end); %% use smallest eigenvalue
    
    cinv=inv(covtrial); %% get inverse
    
    
    cond1_ind=find(trialtypes==1);
    cond2_ind=find(trialtypes==2);
    H_tstat=zeros(length(grid.inside),1);
    tstat=zeros(length(grid.inside),1);
    nx=length(cond1_ind);
    ny=length(cond2_ind);
    dfe = nx + ny - 2;
    power_trial=zeros(Ntrials,length(freq_ind));
    
    for i=1:length(grid.inside),
        lf=cell2mat(grid.leadfield(grid.inside(i)));
        
        %% get optimal orientation- direct copy from Robert's beamformer_lcmv.m
        projpower_vect=pinv(lf'*cinv*lf);
        [u, s, v] = svd(real(projpower_vect));
        eta = u(:,1);
        lf  = lf * eta; %% now have got the lead field at this voxel, compute some contrast
        weights=lf'*cinv/(lf'*cinv*lf); %% Corrected weights calc
        
        if S.return_weights
            ctf_weights(i,:)=weights;
        end
        
        for j=1:Ntrials, %% this non-linear step (power estimation) has to be done at each location
            fdatatrial=squeeze(fftnewdata(j,freq_ind,:))*weights';
            if S.logflag,
                power_trial(j,:)=log(fdatatrial.*conj(fdatatrial));
            else
                power_trial(j,:)=fdatatrial.*conj(fdatatrial); %%
            end; % if
            if S.rankflag,
                power_trial(j,:)=tiedrank(power_trial(j,:));
            end;
            
        end; % for j
        
        xba=mean(power_trial(cond1_ind,:)); %% average power in each freq(across epochs) in cond 1
        yba=mean(power_trial(cond2_ind,:)); %%
        Scov_x=cov(power_trial(cond1_ind,:))*nx; %% covariance of fft bins over epochs
        Scov_y=cov(power_trial(cond2_ind,:))*ny; %% channel covariance
        Scov_pooled=(Scov_x+Scov_y)/(nx+ny-2);
        invS=pinv(Scov_pooled);
        H_tstat(i)=Ntrials*(xba-yba)*invS*(xba-yba)';
        
        %% also compute normal tstat (based on mean power)
        pdiff=mean(xba)-mean(yba);
        xba_epochs=mean(power_trial(cond1_ind,:),2); %% average power (in all frequency) across epochs in cond 1
        yba_epochs=mean(power_trial(cond2_ind,:),2); %%
        s2x=var(xba_epochs);
        s2y=var(yba_epochs);
        sPooled = sqrt(((nx-1) .* s2x + (ny-1) .* s2y) ./ dfe);
        se = sPooled .* sqrt(1./nx + 1./ny);
        tstat(i) = pdiff ./ se; %
        if ~S.logflag,
            normdiff(i)=pdiff./(weights*noise*weights'); %% normalized difference
        else
            normdiff(i)=pdiff; %% if log transform is being used don't need to normalize (as it is a ratio already)
        end; % S
        
        
        if i/100==floor(i/100)
            disp(sprintf('done Hotelling t stats for %3.2f percent of freq band %d of %d, log=%d, rank=%d',100*i/length(grid.inside),f,Nbands,S.logflag,S.rankflag));
        end; % if
    end;
    
    
    %% compute spectra for max voxel for now
    [maxval maxind]=max(H_tstat);
    titlestr='Hotelling';
    [maxval maxind2]=max(tstat);
    titlestr=strvcat(titlestr,'T stat on mean');
    maxinds=[maxind maxind2];
    for i=1:length(maxinds),
        lf=cell2mat(grid.leadfield(grid.inside(maxinds(i))));
        %% get optimal orientation- direct copy from Robert's beamformer_lcmv.m
        projpower_vect=pinv(lf'*cinv*lf);
        [u, s, v] = svd(real(projpower_vect));
        eta = u(:,1);
        lf  = lf * eta; %% now have got the lead field at this voxel, compute some contrast
        weights=(lf'*cinv*lf)*lf'*cinv; %% no regularisation for now
        
        for j=1:Ntrials, %% this non-linear step (power estimation) has to be done at each location
            fdatatrial=squeeze(fftnewdata(j,freq_ind,:))*weights';
            power_trial(j,:)=fdatatrial.*conj(fdatatrial); %% could log transform to make more normal
        end; % for j
        
        xba=mean(power_trial(cond1_ind,:)); %% average change across each condition
        yba=mean(power_trial(cond2_ind,:)); %% average change across each condition
        figure;
        plot(fHz(freq_ind),xba,fHz(freq_ind),yba);
        title(deblank(titlestr(i,:)));
        legend('cond1','cond2');
    end; % for i
    
    
    
    sMRI = fullfile(spm('dir'), 'canonical', 'single_subj_T1.nii');
    
    
    csource=grid;
    csource.pow_H_tstat(csource.inside) = sqrt(H_tstat);
    csource.pow_tstat(csource.inside) = tstat;
    csource.pow_diff(csource.inside)=normdiff;
    csource.pow_diff(csource.outside)=0;
    csource.pow_H_tstat(csource.outside)=0;
    csource.pow_tstat(csource.outside)=0;
    csource.pos = spm_eeg_inv_transform_points(D.inv{D.val}.datareg.toMNI, csource.pos);
    
    % CTF positions inside head
    ctf_inside = spm_eeg_inv_transform_points(D.inv{D.val}.datareg.toMNI, csource.pos(csource.inside,:));
    
    
    
    cfg1 = [];
    cfg1.sourceunits   = 'mm';
    cfg1.parameter = 'pow_H_tstat';
    cfg1.downsample = 1;
    sourceint_H_tstat = ft_sourceinterpolate(cfg1, csource, sMRI);
    
    
    if ~S.rankflag, %% only need to write t images if data hasn't been ranked (ie. on average all the same)
        
        cfg1 = [];
        cfg1.sourceunits   = 'mm';
        cfg1.parameter = 'pow_tstat';
        cfg1.downsample = 1;
        sourceint_tstat = ft_sourceinterpolate(cfg1, csource, sMRI);
        
        cfg1 = [];
        cfg1.sourceunits   = 'mm';
        cfg1.parameter = 'pow_diff';
        cfg1.downsample = 1;
        sourceint_pow_diff = ft_sourceinterpolate(cfg1, csource, sMRI);
    end; % if rankflag
    
    
    %%
    %S.preview=0;
    
    if (isfield(S, 'preview') && S.preview)
        
        cfg1 = [];
        cfg1.funparameter = 'pow_H_tstat';
        cfg1.funcolorlim = [max(csource.pow_H_tstat)/2 max(csource.pow_H_tstat)];
        cfg1.interactive = 'yes';
        figure
        ft_sourceplot(cfg1,sourceint_H_tstat);
        if ~S.rankflag,
            
            cfg1 = [];
            cfg1.funparameter = 'pow_tstat';
            cfg1.funcolorlim = [min(csource.pow_tstat) max(csource.pow_tstat)];
            cfg1.interactive = 'yes';
            figure
            ft_sourceplot(cfg1,sourceint_tstat);
        end;
        
    end
    
    
    %% else %% write out the data sets
    disp('writing images');
    
    cfg = [];
    cfg.sourceunits   = 'mm';
    cfg.parameter = 'pow';
    cfg.downsample = 1;
    % write t stat
    dirname='mvBf_images';
    if S.logflag,
        dirname=[dirname '_log'];
    end;
    if S.rankflag,
        dirname=[dirname '_rank'];
    end;
    
    res = mkdir(D.path, dirname);
    outvol = spm_vol(sMRI);
    outvol.dt(1) = spm_type('float32');
    outvol.fname= fullfile(D.path, dirname, ['Ht_pw_' spm_str_manip(D.fname, 'r') '_' contrast_str '_' num2str(freqbands(f,1)) '-' num2str(freqbands(f,2)) 'Hz' '.nii']);
    outfilenames=strvcat(outvol.fname);
    outvol = spm_create_vol(outvol);
    spm_write_vol(outvol, sourceint_H_tstat.pow_H_tstat);
    if ~S.rankflag,
        %% normal t stat
        outvol = spm_vol(sMRI);
        outvol.dt(1) = spm_type('float32');
        outvol.fname= fullfile(D.path, dirname, ['t_pw_' spm_str_manip(D.fname, 'r') '_' contrast_str '_' num2str(freqbands(f,1)) '-' num2str(freqbands(f,2)) 'Hz' '.nii']);
        outvol = spm_create_vol(outvol);
        outfilenames=strvcat(outvol.fname);
        spm_write_vol(outvol, sourceint_tstat.pow_tstat);
        %% write difference in weight normalised power
        outvol = spm_vol(sMRI);
        outvol.dt(1) = spm_type('float32');
        outvol.fname= fullfile(D.path, dirname, ['N_pw_' spm_str_manip(D.fname, 'r') '_' contrast_str '_' num2str(freqbands(f,1)) '-' num2str(freqbands(f,2)) 'Hz' '.nii']);
        outvol = spm_create_vol(outvol);
        outfilenames=strvcat(outvol.fname);
        spm_write_vol(outvol, sourceint_pow_diff.pow_diff);
    end; % if S.rankflag
    
    
end; % for f
    
end % function


