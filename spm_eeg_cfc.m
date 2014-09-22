function spm_eeg_cfc(S)

% computes GLM for phase-amplitude and amplitude-amplitude coupling as described in van Wijk et al. (submitted)
% Xamp = independent variable to be explained: Xamp = B1*sin(Xphase) + B2*cos(Xphase) + B3*Xlowamp
% additional regressors may be included
% - overall estimates of PAC & AMP are obtained from continuous (or concatenated) data
% - statistical inference of these estimates is performed by dividing the continuous time series into shorter epochs
% - function writes out images of the estimad PAC & AMP, as well as their p values


% filename              = name SPM file containing raw or preprocessed data, either continuous (preferred) or epoched
% channel               = channel name for which the GLM should be calculated
% writeepochs           = write out images for each epoch (0 = no, 1 = yes), default = 0
% window                = duration individual epochs [s] - only needed for continuous data, default = []

% Xamp.Famp             = Centre frequencies for amplitude frequency, e.g. [60:5:100]
% Xamp.Famp_width       = Bandwidth for bandpass filter: Famp +- Famp_width, e.g [10]
% Xamp.Famp_order       = filter order, default = 4

% Rpac.include          = include sin & cos regressors for PAC caluclation (0 = no, 1 = yes)
% Rpac.alpha            = alpha-level for statistical inference, default = 0.05
% Rpac.Fphase           = Centre frequencies for phase frequency, e.g. [5:1:10]
% Rpac.Fphase_width     = Bandwidth for bandpass filter: Fphase +- Fphase_width
% Rpac.Fphase_order     = filter order, default = 2

% Ramp.include          = include low-frequency amplitude regressor for detection amplitude-amplitude coupling (0 = no, 1 = yes)
% Ramp.alpha            = alpha-level for statistical inference, default = 0.05
% Ramp.Flowamp          = Centre frequencies for phase frequency, e.g. [5:1:10]
% Ramp.Flowamp_width    = Bandwidth for bandpass filter: Fphase +- Fphase_width
% Ramp.Flowamp_order    = filter order, default = 2

% custom_reg1.name       = name of additional regressor
% custom_reg1.x          = time series of regressor

% ...

% Rpac.Fphase & Ramp.Flowamp should be equal if both are included in the GLM
% Bernadette van Wijk, Vladimir Litvak

%{
window = 2;
Rpac.include = 1;
Ramp.include = 1;
writeepochs =1;

hamp = 60:10:150;
lamp = 4:2:30;

try Xamp.Famp;              catch, Xamp.Famp = hamp; end
try Xamp.Famp_width;        catch, Xamp.Famp_width = 10; end
try Xamp.Famp_order;        catch, Xamp.Famp_order=4; end
if Rpac.include
    try Rpac.alpha;         catch, Rpac.alpha=.05; end
    try Rpac.Fphase;        catch, Rpac.Fphase = lamp; end
    try Rpac.Fphase_width;  catch, Rpac.Fphase_width=2; end
    try Rpac.Fphase_order;  catch, Rpac.Fphase_order=2; end
end
if Ramp.include
    try Ramp.alpha;         catch, Ramp.alpha=.05; end
    try Ramp.Flowamp;       catch, Ramp.Flowamp = lamp; end
    try Ramp.Flowamp_width; catch, Ramp.Flowamp_width=2; end
    try Ramp.Flowamp_order; catch, Ramp.Flowamp_order=2; end
end
%}

SVNrev = '$Rev: 6186 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','Cross-frequency coupling'); spm('Pointer','Watch');

D = spm_eeg_load(S.D);

if ~isfield(S, 'conditions') || isempty(S.conditions),  S.conditions = D.condlist;  end
if ~iscell(S.conditions), S.conditions = {S.conditions};                            end


if ~isequal(D.transformtype, 'TF')
    error('The input should time-frequency dataset');
end

allamp   = [];
allphase = [];
for i = 1:numel(S.regressors)
    fun  = char(fieldnames(S.regressors{i}));
    S1   = S.regressors{i}.(fun);
    S1.D = D;
    S1.summarise = false;
    res =  feval(['spm_eeg_regressors_' fun], S1);
    
    switch fun
        case 'tfpower'
            allamp   = spm_cat_struct(allamp, res);
        case 'tfphase'
            allphase = spm_cat_struct(allphase, res);
    end
end

allconfounds = [];
for i = 1:numel(S.confounds)
    fun  = char(fieldnames(S.confounds{i}));
    S1   = S.confounds{i}.(fun);
    S1.D = D;
    S1.summarise = false;
    res =  feval(['spm_eeg_regressors_' fun], S1);
    
    allconfounds   = spm_cat_struct(allconfounds, res);
end


freqind = D.indfrequency(min(S.freqwin)):D.indfrequency(max(S.freqwin));
if isempty(freqind) || any(isnan(freqind))
    error('Selected frequency window is invalid.');
end

data = spm_squeeze(mean(D(D.selectchannels(S.channels), freqind, :, D.indtrial(S.conditions, 'GOOD')), 1), 1);
cut  = round(D.fsample/4); %removed at start and end of each filter time series to avoid filter ringing - for trial type data this means a loss of samples per trial


if size(data,3)>1
    datatype     = 'trials';
    trialsamples = size(data, 2)-2*cut+1;
    nepochs      = size(data,3);
    totalsamples = trialsamples*nepochs;
    disp(['number of epochs used for statistics: ', num2str(nepochs)]);
else
    datatype     = 'continuous';
    totalsamples = size(data, 2)-2*cut+1;
    trialsamples = round(window*D.fsample);
    nepochs      = floor(totalsamples/trialsamples);
    disp(['number of epochs used for statistics: ', num2str(nepochs)]);
end

% --------------------------------------------
% get amplitude timeseries
Famp = D.frequencies(freqind);
for N = 1:length(freqind)
    fprintf('\nF amp = %.1f |\t', Famp(N));
    
    if strcmp(datatype,'trials')
        for k = 1:size(data, 3)
            amp_high = data(N, cut:end-cut, k);
            AMP(N,(k-1)*trialsamples+1:k*trialsamples) = amp_high;
            amp(N,k,:) = (amp_high - mean(amp_high))./std(amp_high);
        end
        AMP(N,:) = (AMP(N,:) - mean(AMP(N,:)))./std(AMP(N,:));
        
    elseif strcmp(datatype,'continuous')
        
        amp_high = data(N, cut:end-cut);
        AMP(N,:) = amp_high;
        for k = 1:nepochs
            amp(N,k,:) = AMP(N,(k-1)*trialsamples+1:k*trialsamples);
            amp(N,k,:) = (amp(N,k,:)-mean(amp(N,k,:)))./std(amp(N,k,:));
        end
        AMP(N,:) = (AMP(N,:)-mean(AMP(N,:)))./std(AMP(N,:));
    end
    
end

nsamples = size(data, 2);

% --------------------------------------------
% get phase time series
SINE = {};
sine = {};
COSINE = {};
cosine = {};

for i = 1:numel(allphase)
    PHASE = allphase(i).R;
    
    nphase = 0.5*size(PHASE, 2);
    
    for j = 1:nphase
        
        ind = max(strfind(allphase(i).names{j}, '_'));
        phasefreq(i, j) = sscanf(allphase(1).names{j}(ind+1:end), '%fHz');
        fprintf('\nF phase = %.1f |\t',  phasefreq(i, j));
        
        if strcmp(datatype,'trials')
            for k = 1:size(data, 3)
                phase_low = PHASE(((k-1)*nsamples+1):k*nsamples, j) ;
                SINE{i}(j,(k-1)*trialsamples+1:k*trialsamples) = phase_low(cut:end-cut);
                sine{i}(j,k,:) = phase_low(cut:end-cut);
                
                phase_low = PHASE(((k-1)*nsamples+1):k*nsamples, j + nphase) ;
                COSINE{i}(j,(k-1)*trialsamples+1:k*trialsamples) = phase_low(cut:end-cut);
                cosine{i}(j,k,:) = phase_low(cut:end-cut);
            end
            
        elseif strcmp(datatype,'continuous')
            phase_low = PHASE(:, j) ;
            SINE{i}(j,:) = phase_low(cut:end-cut);
            
            phase_low = PHASE(:, j+nphase) ;
            COSINE{i}(j,:) = phase_low(cut:end-cut);
            for k = 1:nepochs
                sine{i}(j, k,:)   = SINE(j,(k-1)*trialsamples+1:k*trialsamples);
                cosine{i}(j, k,:) = COSINE(j,(k-1)*trialsamples+1:k*trialsamples);
            end
        end
    end
end

% --------------------------------------------
% get amplitude time series for low frequencies
AMP_LOW = {};
amp_low = {};
for i = 1:numel(allamp)
    
    namp = size(allamp(i).R, 2);
    
    for j = 1:namp
        ind = max(strfind(allamp(i).names{j}, '_'));
        ampfreq(i, j) =  sscanf(allamp(1).names{j}(ind+1:end), '%fHz');
        fprintf('\nF low amp = %.1f |\t', ampfreq(i, j));
        
        if strcmp(datatype,'trials')
            for k = 1:size(data, 3)
                
                amplow = allamp(i).R(((k-1)*nsamples+1):k*nsamples, j);
                AMP_LOW{i}(j,(k-1)*trialsamples+1:k*trialsamples)=amplow(cut:end-cut);
                amp_low{i}(j,k,:)=(amplow(cut:end-cut)-mean(amplow(cut:end-cut)))./std(amplow(cut:end-cut));
                
            end
            AMP_LOW{i}(j,:)=(AMP_LOW{i}(j,:)-mean(AMP_LOW{i}(j,:)))./std(AMP_LOW{i}(j,:));
            
        elseif strcmp(datatype,'continuous')
            
            amplow= allamp(i).R(:, j);
            AMP_LOW{i}(j,:)=amplow(cut:end-cut);
            for k=1:nepochs
                amp_low{i}(j,k,:)=AMP_LOW{i}(j,(k-1)*trialsamples+1:k*trialsamples);
                amp_low{i}(j,k,:)=(amp_low{i}(j,k,:)-mean(amp_low{i}(j,k,:)))./std(amp_low{i}(j,k,:));
            end
            AMP_LOW{i}(j,:)=(AMP_LOW{i}(j,:)-mean(AMP_LOW{i}(j,:)))./std(AMP_LOW{i}(j,:));
        end
    end
end


% --------------------------------------------
% get amplitude time series for low frequencies
CONFOUNDS = {};
confounds = {};
for i = 1:numel(allconfounds)
    
    nconf = size(allconfounds(i).R, 2);
    
    for j = 1:nconf
        fprintf('\nConfound: %s |\t', allconfounds(i).names{j});
        
        if strcmp(datatype,'trials')
            for k = 1:size(data, 3)
                
                conf = allconfounds(i).R(((k-1)*nsamples+1):k*nsamples, j);
                CONFOUNDS{i}(j,(k-1)*trialsamples+1:k*trialsamples) = conf(cut:end-cut);
                confounds{i}(j,k,:) = (conf(cut:end-cut)-mean(conf(cut:end-cut)))./std(conf(cut:end-cut));
                
            end
            CONFOUNDS{i}(j,:)=(CONFOUNDS{i}(j,:)-mean(CONFOUNDS{i}(j,:)))./std(CONFOUNDS{i}(j,:));
            
        elseif strcmp(datatype,'continuous')
            
            conf = allconfounds(i).R(:, j);
            CONFOUNDS{i}(j,:) = conf(cut:end-cut);
            for k=1:nepochs
                confounds{i}(j,k,:) = CONFOUNDS{i}(j,(k-1)*trialsamples+1:k*trialsamples);
                confounds{i}(j,k,:)=(confounds{i}(j,k,:)-mean(confounds{i}(j,k,:)))./std(confounds{i}(j,k,:));
            end
            CONFOUNDS{i}(j,:) = (CONFOUNDS{i}(j,:)-mean(CONFOUNDS{i}(j,:)))./std(CONFOUNDS{i}(j,:));
        end
    end
end

fprintf('\n\n')

Flow = cat(1, phasefreq, ampfreq);
% Could this possibly be relaxed?
if any(any(diff(Flow, [], 1)))
    error('The frequency axes for all regressors should be identical');
else
    Flow = Flow(1, :);
end

% --------------------------------------------
% compute GLM

for j=1:length(Flow)
    fprintf('%d  ',Flow(j))
    for N=1:length(Xamp.Famp)
        
        % --------------------------------------------
        % GLM for all data appended
        
        if Rpac.include==1 && Ramp.include == 1
            
            % [additional custom regressors to be added]
            X=[(cos(PHASE(j,:))-mean(cos(PHASE(j,:))))./std(cos(PHASE(j,:)));(sin(PHASE(j,:))-mean(sin(PHASE(j,:))))./std(sin(PHASE(j,:)));AMP_LOW(j,:)]';
            y=AMP(N,:);
            
            V=[];
            c=[1;1;1];
            [T,df,all_Beta(N,j,:),xX,xCon]=spm_ancova(X,V,y',c);
            
            all_SSy(N,j)=sum((y-mean(y)).^2);
            all_residuals=y-(all_Beta(N,j,1).*X(:,1)+all_Beta(N,j,2).*X(:,2))';
            all_SSe(N,j)=sum((all_residuals-mean(all_residuals)).^2);
            all_r(N,j)=real(sqrt((all_SSy(N,j)-all_SSe(N,j))/all_SSy(N,j)));
            
            all_residuals_amp=y-(all_Beta(N,j,3).*X(:,3))';
            all_SSe_amp(N,j)=sum((all_residuals_amp-mean(all_residuals_amp)).^2);
            all_r_amp(N,j)=real(sqrt((all_SSy(N,j)-all_SSe_amp(N,j))/all_SSy(N,j)));
            
            all_residuals_total=y-(all_Beta(N,j,1).*X(:,1)+all_Beta(N,j,2).*X(:,2)+all_Beta(N,j,3).*X(:,3))';
            all_SSe_total(N,j)=sum((all_residuals_total-mean(all_residuals_total)).^2);
            all_r_total(N,j)=real(sqrt((all_SSy(N,j)-all_SSe_total(N,j))/all_SSy(N,j)));
            
            all_Bpac(N,j)=sqrt(all_Beta(N,j,1).^2+all_Beta(N,j,2).^2);
            all_Bamp(N,j)=all_Beta(N,j,3);
            
        elseif Rpac.include==1 && Ramp.include == 0
            
            % [additional custom regressors to be added]
            X=[(cos(PHASE(j,:))-mean(cos(PHASE(j,:))))./std(cos(PHASE(j,:)));(sin(PHASE(j,:))-mean(sin(PHASE(j,:))))./std(sin(PHASE(j,:)))]';
            y=AMP(N,:);
            
            V=[];
            c=[1;1];
            [T,df,all_Beta(N,j,:),xX,xCon]=spm_ancova(X,V,y',c);
            
            all_SSy(N,j)=sum((y-mean(y)).^2);
            all_residuals=y-(all_Beta(N,j,1).*X(:,1)+all_Beta(N,j,2).*X(:,2))';
            all_SSe(N,j)=sum((all_residuals-mean(all_residuals)).^2);
            all_r(N,j)=real(sqrt((all_SSy(N,j)-all_SSe(N,j))/all_SSy(N,j)));
            
            all_Bpac(N,j)=sqrt(all_Beta(N,j,1).^2+all_Beta(N,j,2).^2);
            
        elseif Rpac.include==0 && Ramp.include == 1
            
            % [additional custom regressors to be added]
            X=[AMP_LOW(j,:)]';
            y=AMP(N,:);
            
            V=[];
            c=[1];
            [T,df,all_Beta(N,j,:),xX,xCon]=spm_ancova(X,V,y',c);
            
            all_SSy(N,j)=sum((y-mean(y)).^2);
            all_residuals_amp=y-(all_Beta(N,j,1).*X(:,1))';
            all_SSe_amp(N,j)=sum((all_residuals_amp-mean(all_residuals_amp)).^2);
            all_r_amp(N,j)=real(sqrt((all_SSy(N,j)-all_SSe_amp(N,j))/all_SSy(N,j)));
            
            all_Bamp(N,j)=all_Beta(N,j,1);
            
        end
        
        % --------------------------------------------
        % GLM per trial
        
        for k=1:size(data,2)
            
            if Rpac.include==1 && Ramp.include == 1
                Xk=[(cos(squeeze(phase(j,k,:)))-mean(cos(squeeze(phase(j,k,:)))))./std(cos(squeeze(phase(j,k,:)))),(sin(squeeze(phase(j,k,:)))-mean(sin(squeeze(phase(j,k,:)))))./std(sin(squeeze(phase(j,k,:)))),squeeze(amp_low(j,k,:))]';
                yk=squeeze(amp(N,k,:));
                [T,df,Beta(:,k),xX,xCon]=spm_ancova(Xk',V,yk,c); %to get betas
                if writeepochs
                    Beta1(N,j,k)=Beta(1,k);
                    Beta2(N,j,k)=Beta(2,k);
                    Beta3(N,j,k)=Beta(3,k);
                end
            elseif Rpac.include==1 && Ramp.include == 0
                Xk=[(cos(squeeze(phase(j,k,:)))-mean(cos(squeeze(phase(j,k,:)))))./std(cos(squeeze(phase(j,k,:)))),(sin(squeeze(phase(j,k,:)))-mean(sin(squeeze(phase(j,k,:)))))./std(sin(squeeze(phase(j,k,:))))]';
                yk=squeeze(amp(N,k,:));
                [T,df,Beta(:,k),xX,xCon]=spm_ancova(Xk',V,yk,c); %to get betas
                if writeepochs
                    Beta1(N,j,k)=Beta(1,k);
                    Beta2(N,j,k)=Beta(2,k);
                end
            elseif Rpac.include==0 && Ramp.include == 1
                Xk=[squeeze(amp_low(j,k,:))]';
                yk=squeeze(amp(N,k,:));
                [T,df,Beta(3,k),xX,xCon]=spm_ancova(Xk',V,yk,c); %to get betas
                if writeepochs
                    Beta3(N,j,k)=Beta(3,k);
                end
            end
        end %k
        
        % --------------------------------------------
        % test for significance
        
        if Rpac.include
            Xb=[];
            Xb(1:k,1)=ones(k,1);Xb(k+1:2*k,2)=ones(k,1);
            yb=[Beta(1,:),Beta(2,:)];
            c=[1;1];
            [Tb,df,Beta_b,xX,xCon]=spm_ancova(Xb,V,yb',c);
            F=Tb^2;
            pb(N,j)=1-spm_Fcdf(F,df(1),df(2));
        end
        
        if Ramp.include
            yb=[Beta(3,:)];
            % [H,P] = ttest(yb);
            % pb_amp(N,j)=P;
        end
        
        if Rpac.include && Ramp.include
            Xb=[];
            Xb(1:k,1)=ones(k,1);Xb(k+1:2*k,2)=ones(k,1);Xb(2*k+1:3*k,3)=ones(k,1);
            c=[1;1;1];
            yb=[Beta(1,:),Beta(2,:),Beta(3,:)];
            [Tb,df,Beta_b,xX,xCon]=spm_ancova(Xb,V,yb',c);
            F_total=Tb^2;
            pb_total(N,j)=1-spm_Fcdf(F_total,df(1),df(2));
        end
        
    end %N
end %j

sig_r=(pb<=Rpac.alpha);
sig_r_amp=0;%(pb_amp<=Ramp.alpha);
sig_r_total=(pb_total<=Rpac.alpha);

figure,
if Rpac.include
    subplot(3,2,1),imagesc(Flow,Xamp.Famp,all_r),set(gca,'ydir','normal');title(['r PAC']);colorbar;
    subplot(3,2,2),imagesc(Flow,Xamp.Famp,sig_r),set(gca,'ydir','normal');title(['significant r PAC']), colorbar;
end
if Ramp.include
    subplot(3,2,3),imagesc(Flow,Xamp.Famp,all_r_amp),set(gca,'ydir','normal');title(['c AMP']);colorbar;
    subplot(3,2,4),imagesc(Flow,Xamp.Famp,sig_r_amp),set(gca,'ydir','normal');title(['significant c AMP']), colorbar;
end
if Rpac.include && Ramp.include
    subplot(3,2,5),imagesc(Flow,Xamp.Famp,all_r_total),set(gca,'ydir','normal');title(['r total']);colorbar;
    subplot(3,2,6),imagesc(Flow,Xamp.Famp,sig_r_total),set(gca,'ydir','normal');title(['significant r total']), colorbar;
end

% --------------------------------------------
% write out images

cnt=1;

if Rpac.include
    image(cnt).val     = all_r;
    image(cnt).label   = ['r_pac_'  spm_file(fname, 'basename')];
    cnt=cnt+1;
    image(cnt).val     = pb;
    image(cnt).label   = ['p_pac_'  spm_file(fname, 'basename')];
    cnt=cnt+1;
    image(cnt).val     = sig_r;
    image(cnt).label   = ['sig_pac_'  spm_file(fname, 'basename')];
    cnt=cnt+1;
    image(cnt).val     = squeeze(all_Beta(N,j,1));
    image(cnt).label   = ['r_B1_'  spm_file(fname, 'basename')];
    cnt=cnt+1;
    image(cnt).val     = squeeze(all_Beta(N,j,2));
    image(cnt).label   = ['r_B2_'  spm_file(fname, 'basename')];
    cnt=cnt+1;
end

if Ramp.include
    % image(cnt).val     = all_r_amp;
    % image(cnt).label   = ['r_amp_'  spm_file(fname, 'basename')];
    % cnt=cnt+1;
    image(cnt).val     = all_Bamp;
    image(cnt).label   = ['c_amp_'  spm_file(fname, 'basename')];
    cnt=cnt+1;
    image(cnt).val     = pb_amp;
    image(cnt).label   = ['p_amp_'  spm_file(fname, 'basename')];
    cnt=cnt+1;
    image(cnt).val     = sig_r_amp;
    image(cnt).label   = ['sig_amp_'  spm_file(fname, 'basename')];
    cnt=cnt+1;
end

if Rpac.include && Ramp.include
    image(cnt).val     = all_r_total;
    image(cnt).label   = ['r_total_'  spm_file(fname, 'basename')];
    cnt=cnt+1;
    image(cnt).val     = pb_total;
    image(cnt).label   = ['p_total_'  spm_file(fname, 'basename')];
    cnt=cnt+1;
    image(cnt).val     = sig_r_total;
    image(cnt).label   = ['sig_total_'  spm_file(fname, 'basename')];
    cnt=cnt+1;
end

%per trial
if writeepochs
    for k=1:ntrials
        if Rpac.include
            image(cnt).val = squeeze(Beta1(:,:,k));
            image(cnt).label   = ['trial',num2str(k),'_B1_'  spm_file(fname, 'basename')];
            cnt=cnt+1;
            image(cnt).val = squeeze(Beta2(:,:,k));
            image(cnt).label   = ['trial',num2str(k),'_B2_'  spm_file(fname, 'basename')];
            cnt=cnt+1;
        end
        if Ramp.include
            image(cnt).val = squeeze(Beta3(:,:,k));
            image(cnt).label   = ['trial',num2str(k),'_B3_'  spm_file(fname, 'basename')];
            cnt=cnt+1;
        end
    end
end


%% WRITE IMAGES

% ...
% ...
% ...

