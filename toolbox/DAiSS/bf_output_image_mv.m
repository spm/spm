function res = bf_output_image_mv(BF, S)
% Computes multivariate test on a number of frequency bands
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes, modified from Vladimir Litvak's example code
% $Id: bf_output_image_mv.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------

%% no covariance matrix creation -so need to check bands are within this window

if nargin == 0
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
    
    %     design = cfg_const;
    %     design.tag = 'design';
    %     design.name = 'design';
    %     design.help = {'Use default settings for the inversion'};
    %     design.val  = {1};
    
    design = cfg_files;
    design.tag = 'design';
    design.name = 'design matrix';
    design.filter = 'mat';
    design.num = [1 Inf];
    design.help = {'Select the design matrix'};
    
    
    
    woi = cfg_entry;
    woi.tag = 'woi';
    woi.name = 'Time windows of interest';
    woi.strtype = 'r';
    woi.num = [Inf 2];
    woi.help = {'Time windows (in ms)'};
    
    foi = cfg_entry;
    foi.tag = 'foi';
    foi.name = 'Frequency bands of interest';
    foi.strtype = 'r';
    foi.num = [Inf 2];
    foi.help = {'Freq bands (in Hz)'};
    
    
    datafeatures = cfg_menu;
    datafeatures.tag = 'datafeatures';
    datafeatures.name = 'Features';
    datafeatures.labels =get_data_features;
    datafeatures.values =get_data_features;
    
    datafeatures.help = {'Data features of interest'};
    datafeatures.val =datafeatures.values(end);
    
    contrast = cfg_entry;
    contrast.tag = 'contrast';
    contrast.name = 'contrast';
    contrast.strtype = 'i';
    contrast.num = [1 Inf];
    contrast.val = {1};
    
    result         = cfg_menu;
    result.tag     = 'result';
    result.name    = 'What to output';
    result.help    = {'Specify output type'};
    result.labels  = {
        'chi square'
        'BIC'
        'r square'
        }';
    result.values  = result.labels;
    result.val = {'chi square'};
    
    modality         = cfg_menu;
    modality.tag     = 'modality';
    modality.name    = 'Modality';
    modality.help    = {'Specify modality'};
    modality.labels  = {
        'MEG'
        'MEGPLANAR'
        'EEG'
        }';
    modality.values  = {
        'MEG'
        'MEGPLANAR'
        'EEG'
        }';
    modality.val = {'MEG'};
    
    custom = cfg_branch;
    custom.tag = 'custom';
    custom.name = 'Custom';
    custom.help = {'Define custom settings for the inversion'};
    custom.val  = {whatconditions, contrast, woi};
    
    isdesign = cfg_choice;
    isdesign.tag = 'isdesign';
    isdesign.name = 'Design matrix parameters';
    isdesign.help = {'Choose whether to load custom design'};
    isdesign.values = {design, custom};
    isdesign.val = {design};
    
    sametrials = cfg_menu;
    sametrials.tag = 'sametrials';
    sametrials.name = 'Trials same as for filters';
    sametrials.labels = {'yes', 'no'};
    sametrials.values = {true, false};
    sametrials.val = {false};
    sametrials.help = {'Take the same trials as used for filter computation',...
        'This is useful for bootstrap.'};
    
    image_mv      = cfg_branch;
    image_mv.tag  = 'image_mv';
    image_mv.name = 'Mv image';
    image_mv.val  = { isdesign, datafeatures, foi,  result, sametrials, modality};
    
    res = image_mv;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end





D = BF.data.D;


if isfield(S.isdesign,'custom'),
    %% gui specified conditions and contrast
        
         
    woitmp = S.isdesign.custom.woi;
   
    % DP - The indsample function does not work for matrices, so I have
    % looped through. Otherwise, we are left with only one windows of
    % interest.
    for wi = 1:size(woitmp,1)
        woiind(wi,:)=D.indsample(woitmp(wi,:)/1000);
        woi(wi,:)=D.time(woiind(wi,:)); %% in seconds
    end
   
    duration=unique(woiind(:,2)-woiind(:,1))./D.fsample; %% in sec
    if numel(duration)>1,
        error('both windows need to be the same length');
    end;
   
    % DP - this can result in a tiny bit of residual, presumably due to
    % numerical imprecision. Not sure. Doesn't seem necessary anyway since
    % the number of samples is equal.
%     duration=unique(woi(:,2)-woi(:,1));
%     if numel(duration)>1,
%         error('both windows need to be the same length');
%     end;
%     woi = S.isdesign.custom.woi;
%     woiind=D.indsample(woi/1000);
%     woi=D.time(woiind); %% in seconds
    
%     duration=unique(woiind(:,2)-woiind(:,1))./D.fsample; %% in sec
%     if numel(duration)>1,
%         error('both windows need to be the same length');
%     end;
    
    duration=unique(woi(:,2)-woi(:,1));
    if numel(duration)>1,
        error('both windows need to be the same length');
    end;
    
    %
    whatconditions=S.isdesign.custom.whatconditions;
    
    if isfield(whatconditions, 'all')
        if S.sametrials
            trials{1} = BF.features.trials;
        else
            trials{1} = 1:D.ntrials;
        end
        clabel{1}='all';
    else
        for i = 1:numel(whatconditions.condlabel)
            if isempty(D.indtrial(whatconditions.condlabel{i}, 'GOOD'))
                error('No trials matched the selection.');
            end
            
            clabel{i}=whatconditions.condlabel{i};
            
            if S.sametrials
                trials{i} = BF.features.trials(strmatch(clabel{i}, D.conditions(BF.features.trials)));
            else
                trials{i} = D.indtrial(whatconditions.condlabel{i}, 'GOOD');
            end
                        
            
        end
        if isempty(trials)
            error('No trials matched the selection, check the specified condition labels');
        end
    end
    
      
    %%check for number of trials //// added by ANNA 30/04/2013
    for i=1:numel(trials)
        num_trials(i) = length(trials{i});
    end
    
    if min(num_trials)~= max(num_trials)
        warning ('Number of trials are not the same accross conditions- throwing away');
        num_trials
    end
    nt = min(num_trials); % ANNA -throw away the extra trials
    
    col=0;
    X=[];
    Xtrials=[];
    Xstartlatencies=[];
    xlabel=[];
    for j=1:size(woi, 1),
        for i=1:numel(trials)
            col=col+1;
            %nt=numel(trials{i}); % ANNA look above for the warning
            Xtmp=[zeros(size(X,1),1); ones(nt,1)];
            if col>1,
                X=[X;zeros(nt,size(X,2))]; %1)]; ANNA
            end;
            X=[X Xtmp];
            tlist = spm_vec(trials{i}); % ANNA list of trials 
            Xtrials=[Xtrials ; tlist(1:nt)]; % ANNA list of trials - this must be a vector
            Xstartlatencies=[Xstartlatencies; ones(nt,1).*woi(j,1)];%% again this has to be a vector
            xlabel=strvcat(xlabel,[clabel{i} ',' num2str(woi(j,1))]);
        end;
    end;
%     % ORIG
%     col=0;
%     X=[];
%     Xtrials=[];
%     Xstartlatencies=[];
%     xlabel=[];
%     for j=1:size(woi, 1),
%         for i=1:numel(trials)
%             col=col+1;
%             nt=numel(trials{i});
%             Xtmp=[zeros(size(X,1),1); ones(nt,1)];
%             if col>1,
%                 X=[X;zeros(nt,1)];
%             end;
%             X=[X Xtmp];
%             Xtrials=[Xtrials ;[trials{i}]];
%             Xstartlatencies=[Xstartlatencies; ones(nt,1).*woi(j,1)];
%             xlabel=strvcat(xlabel,[clabel{i} ',' num2str(woi(j,1))]);
%         end;
%     end;
    
    allsamples=[D.indsample(Xstartlatencies); D.indsample(Xstartlatencies+ones(size(Xstartlatencies)).*duration)]';
    contrast=S.isdesign.custom.contrast';
    
    
    
    
    
else %%  conditions and contrast  specified in a file
    
    if ~exist(cell2mat(S.isdesign.design)),
        error('Cannot load design matrix');
    end;
    disp('loading design matrix');
    a=load(cell2mat(S.isdesign.design));
    X=a.design.X; %% design matrix
    contrast=a.design.contrast;
    ntrials=size(X,1);
    if (size(a.design.Xstartlatencies,1)~=ntrials)||(size(a.design.Xtrials,1)~=ntrials)
        error('start latencies and Xtrials and X should have a value per row of the design');
    end;
    
    Xtrials=a.design.Xtrials; %% indices of trials to use
    for j=1:ntrials,
        allsamples(j,1)=D.indsample(a.design.Xstartlatencies(j));
        allsamples(j,2)=D.indsample(a.design.Xstartlatencies(j)+a.design.Xwindowduration);
    end;
    
end;


if size(Xtrials)~=size(Xstartlatencies)
    error('Xtrials and start latencies must be the same length');
end;
if size(Xtrials,1)~=size(X,1)
    error('X must have same number of rows as trials and start latencies');
end;
if size(contrast,1)~=size(X,2),
    error('contrast needs to match number of columns in design');
end;


Fgraph  = spm_figure('GetWin','Graphics'); spm_figure('Clear',Fgraph);

subplot(4,1,1);
imagesc(X);
title('X');
set(gca,'Xtick',1:col);
set(gca,'Xticklabel',xlabel);
subplot(4,1,2);

imagesc(contrast);
title('c');
subplot(4,1,3);

imagesc(Xtrials);
title('trials');
subplot(4,1,4);

imagesc(Xstartlatencies);
title('start latency');

chanind = BF.features.(S.modality).chanind;
U       = BF.features.(S.modality).U;

Nchans=size(U,2); %% effective number of channels

nsamples = unique(allsamples(:,2)-allsamples(:,1));

if length(nsamples) > 1
    error('All time windows should be equal lentgh')
end

alltrials = spm_vec(Xtrials);
ntrials   = length(alltrials);


%% now identify frequency bands of interest

Nbands=size(S.foi,1);
windowduration=(nsamples/D.fsample);
dctfreq = (0:nsamples-1)/2/windowduration;           % DCT frequencies (Hz)
dctT      = spm_dctmtx(nsamples,nsamples);

freqstr=[];
allfreqind=[];
for fband=1:Nbands, %% allows one to break up spectrum and ignore some frequencies
    freqrange=S.foi(fband,:);
    j      = find( (dctfreq >= freqrange(1)) & (dctfreq<=freqrange(2)) );
    featureind{fband}=j;
    allfreqind=sort(unique([allfreqind j]));
    freqstr=[freqstr sprintf('%3.1f-%3.1f,',dctfreq(min(j)),dctfreq(max(j)))];
end; % for fband=1:Nbands

if isempty(allfreqind),
    error('No valid frequency range found');
end;
% Tfull      = dctT(:,allfreqind); %% A filter for all bands (not necessarily continuous)
%% end of freq band section

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', ntrials , 'Computing covariance'); drawnow;
if ntrials  > 100, Ibar = floor(linspace(1, ntrials ,100));
else Ibar = 1:ntrials; end

%% load in data and make up simple design matrix
%ncond=numel(samples); %% number of conditions= columns in design matrix
%Nt=ntrials*ncond; %% total number of time windows under consideration
%X=zeros(Nt,ncond);

flatdata=zeros(ntrials*nsamples,Nchans);
%% want flatdata in form Nchans,Nt*Nsamples

%count=0;
for i = 1:ntrials
    %for j = 1:numel(samples)
    %   count=count+1;
    %   X(count,j)=1;
    Y  = U'*squeeze(D(chanind, allsamples(i,1):allsamples(i,2)-1, alltrials(i)));
    Y  = detrend(Y'); %% detrend and throw away low freq drift
    flatdata((i-1)*nsamples+1:i*nsamples,:) =Y;
    
    %end
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');

W = BF.inverse.(S.modality).W;
nvert = numel(W);
S.regressout=[]; %% turn off for now
regressout=S.regressout;


%% set up the data features
weights=-1; %% set up flag
Yfull=get_data_features(flatdata,nsamples,ntrials,weights,dctT,S.datafeatures,featureind,regressout); %% set up data structures



spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', nvert, 'Scanning grid points'); drawnow;
if nvert > 100, Ibar = floor(linspace(1, nvert,100));
else Ibar = 1:nvert; end

outval = nan(1, nvert);


for i = 1:nvert
    
    if ~isnan(W{i})
        
        w    = W{i};
        
        %% returns columns of a matrix with rows as different observations
        [Yfull,vedata]=get_data_features(flatdata,nsamples,ntrials,w,dctT,S.datafeatures,featureind,regressout); %% extract the data features
        
        Yfull=Yfull-repmat(mean(Yfull),size(Yfull,1),1); %% remove dc level from each column/feature
        Yfull=Yfull./repmat(std(Yfull),size(Yfull,1),1); %% normalize features to have unit variance by default
        [chival,BIC,cva] = output_image_mv_cva(X,Yfull,contrast); %% run the multivariate test
        
        
        switch S.result
            case 'chi square'
                resultstr='chisq';
                outval(i) = chival(1);
            case 'r square'
                outval(i) = cva.ccorr.^2;
                resultstr='rsq';
            case 'BIC'
                outval(i)=BIC(1);
                resultstr='BIC';
        end;
        
    end
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

spm_progress_bar('Clear');


image(1).val   = outval;

image(1).label = ['mv' resultstr S.datafeatures freqstr  spm_file(D.fname, 'basename')];


res = image;