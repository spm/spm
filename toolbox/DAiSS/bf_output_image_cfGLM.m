function res = bf_output_image_cfGLM(BF, S)
% Computes phase-amplitude coupling using a general linear model
% currently takes both low frequency phase and amplitude as regressors
% needs epoched data - uses epochs for statistics
% writes out images for summary phase-amplitude coupling and
% amplitude-amplitude coupling, as well as B coefficients per trial

% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Bernadette van Wijk, Vladimir Litvak
% $Id: bf_output_image_cfGLM.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
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
    
    sametrials = cfg_menu;
    sametrials.tag = 'sametrials';
    sametrials.name = 'Trials same as for filters';
    sametrials.labels = {'yes', 'no'};
    sametrials.values = {true, false};
    sametrials.val = {false};
    sametrials.help = {'Take the same trials as used for filter computation',...
        'This is useful for bootstrap.'};
    
    woi = cfg_entry;
    woi.tag = 'woi';
    woi.name = 'Time window of interest';
    woi.strtype = 'r';
    woi.num = [1 2];
    woi.val = {[-Inf Inf]};
    woi.help = {'Time windows (in ms)'};
    
    phasefreq         = cfg_entry;
    phasefreq.tag     = 'phasefreq';
    phasefreq.name    = 'Phase frequencies';
    phasefreq.strtype = 'r';
    phasefreq.num     = [1 Inf];
    phasefreq.val     = {5:3:30};
    phasefreq.help    = {'Frequencies to compute phase for (as a vector)'};
    
    phaseres         = cfg_entry;
    phaseres.tag     = 'phaseres';
    phaseres.name    = 'Phase resolution';
    phaseres.strtype = 'r';
    phaseres.num     = [1 1];
    phaseres.val     = {2};
    phaseres.help    = {'Frequency resolution for phase computation'};
    
    lowampres         = cfg_entry;
    lowampres.tag     = 'lowampres';
    lowampres.name    = 'Low Amplitude resolution';
    lowampres.strtype = 'r';
    lowampres.num     = [1 1];
    lowampres.val     = {4};
    lowampres.help    = {'Frequency resolution for amplitude computation at low frequencies'};
    
    ampfreq         = cfg_entry;
    ampfreq.tag     = 'ampfreq';
    ampfreq.name    = 'Amplitude frequencies';
    ampfreq.strtype = 'r';
    ampfreq.num     = [1 Inf];
    ampfreq.val     = {30:5:100};
    ampfreq.help    = {'Frequencies to compute amplitude for (as a vector)'};
    
    ampres         = cfg_entry;
    ampres.tag     = 'ampres';
    ampres.name    = 'Amplitude resolution';
    ampres.strtype = 'r';
    ampres.num     = [1 1];
    ampres.val     = {15};
    ampres.help    = {'Frequency resolution for amplitude computation'};
    
    name = cfg_entry;
    name.tag = 'name';
    name.name = 'Channel name';
    name.strtype = 's';
    name.num = [1 Inf];
    name.help = {'Reference channel name.'};
    
    feature         = cfg_menu;
    feature.tag     = 'feature';
    feature.name    = 'Reference feature';
    feature.help    = {'What to take from the reference'};
    feature.labels  = {'Amplitude', 'Phase'};
    feature.values  = {'amplitude', 'phase'};
    
    refchan      = cfg_branch;
    refchan.tag  = 'refchan';
    refchan.name = 'Reference channel';
    refchan.val  = {name, feature};
    
    within = cfg_const;
    within.tag = 'within';
    within.name = 'Within source.';
    within.val  = {1};
    within.help = {'Within source PAC (no reference)'};
    
    reference = cfg_choice;
    reference.tag = 'reference';
    reference.name = 'Reference type';
    reference.values = {within, refchan};
    reference.val = {within};
    
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
    
    ncomp_amp         = cfg_entry;
    ncomp_amp.tag     = 'ncomp_amp';
    ncomp_amp.name    = 'Number of dipole orientations for AMPLITUDE';
    ncomp_amp.strtype = 'r';
    ncomp_amp.num     = [1 1];
    ncomp_amp.val     = {1};
    ncomp_amp.help    = {'Number of dipole orientations for each MEG source for amplitude'};
    
    ncomp_phase         = cfg_entry;
    ncomp_phase.tag     = 'ncomp_phase';
    ncomp_phase.name    = 'Number of dipole orientations for PHASE';
    ncomp_phase.strtype = 'r';
    ncomp_phase.num     = [1 1];
    ncomp_phase.val     = {1};
    ncomp_phase.help    = {'Number of dipole orientations for each MEG source for phase'};
    
    outputname = cfg_entry;
    outputname.tag = 'outputname';
    outputname.name = 'Name output images';
    outputname.strtype = 's';
    outputname.num = [1 Inf];
    outputname.val = {['_']};
    outputname.help = {'To specify details that will be added to the output images file names.'};
    
    image_cfGLM      = cfg_branch;
    image_cfGLM.tag  = 'image_cfGLM';
    image_cfGLM.name = 'cross-frequency GLM image';
    image_cfGLM.val  = {whatconditions, sametrials, woi, phasefreq, ....
        phaseres, lowampres, ampfreq, ampres, reference, modality, ncomp_amp, ncomp_phase, outputname};
    
    res = image_cfGLM;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

ncomponents=S.ncomp_amp*S.ncomp_phase;

D = BF.data.D;

S.woi = 1e-3*S.woi; % ms -> s

samples =  D.indsample(S.woi(1)):D.indsample(S.woi(2));
nsamples = length(samples);
times = D.time(samples);

if isfield(S.whatconditions, 'all')
    S.whatconditions.condlabel = D.condlist;
end

for i = 1:numel(S.whatconditions.condlabel)
    if S.sametrials
        trials{i} = BF.features.trials(strmatch(S.whatconditions.condlabel{i},...
            D.conditions(BF.features.trials)));
    else
        trials{i} = D.indtrial(S.whatconditions.condlabel{i}, 'GOOD');
    end
    
    if isempty(trials{i})
        error('No trials matched the selection.');
    end
    
end

if isempty(trials)
    error('No trials matched the selection, check the specified condition labels');
end

if S.ncomp_amp>1 || S.ncomp_phase>1;
    error('Not supported yet for more than 1 dipole moment.');
end

channels = BF.features.(S.modality).chanind;
U        = BF.features.(S.modality).U;
nchan    = size(U, 2);

alltrials = spm_vec(trials);
ntrials   = length(alltrials);

nphase  = length(S.phasefreq);
namp    = length(S.ampfreq);

W = BF.inverse.(S.modality).W;
nvert = numel(W);

Y = U'*reshape(D(channels, samples, alltrials), nchan, []);
Y = reshape(Y, size(Y, 1), nsamples, ntrials);

if isequal(char(fieldnames(S.reference)), 'refchan')
    ref_feature = S.reference.refchan.feature;
    Yr = D(D.indchannel(S.reference.refchan.name), samples, alltrials);
else
    ref_feature = 'amplitude';
    Yr = Y;
end

switch ref_feature
    case 'amplitude'
        freqoi = S.ampfreq;
        width  = S.ampres;
    case 'phase'
        freqoi = S.phasefreq;
        width  = S.phaseres;
        lowamp_width = S.lowampres;
end

refsig = cell(1, length(freqoi));
for j = 1:length(freqoi)
    refsig{j} = zeros(size(Yr, 1), nsamples-(2*round(D.fsample/6))+1, ntrials);
end

spm_progress_bar('Init', ntrials, ...
    sprintf('Computing reference spectra')); drawnow;
if ntrials > 100, Ibar = floor(linspace(1, ntrials,100));
else Ibar = 1:ntrials; end


for i = 1:ntrials
    spectrum = ft_specest_hilbert(squeeze(Yr(:,:, i)), times,...
        'freqoi', freqoi, 'width', width, 'filttype', 'but', 'filtorder', 2,...
        'filtdir', 'twopass', 'verbose', 0);
    
    for j = 1:length(freqoi)
        tmp = spm_squeeze(spectrum(:, j, :), 2);
        refsig{j}(:,:,i) = tmp(:,round(D.fsample/6):end-round(D.fsample/6)); %to remove edge artefacts
    end
    
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end
spm_progress_bar('Clear');

for j = 1:numel(refsig)
    refsig{j} = reshape(refsig{j}, size(refsig{j}, 1), []);
end

%%%%%%%%%%%%

try S.reference.refchan;
    
    if strmatch(ref_feature,'amplitude');
        
        for j = 1:size(refsig{1},1)%1:numel(refsig)
            amplitude(j,:)=abs(refsig{j});
            AMPLITUDE(j,:)=((amplitude(j,:)-mean(amplitude(j,:)))./std(amplitude(j,:)))';
        end
        
        amplitude=reshape(amplitude,numel(refsig),nsamples-(2*round(D.fsample/6))+1, ntrials);
        
        for j = 1:size(refsig{1},1)% 1:numel(refsig)
            for k= 1:ntrials
                amplitude(j,:,k)=((amplitude(j,:,k)-mean(amplitude(j,:,k)))./std(amplitude(j,:,k)))';
            end
        end
        
        
    elseif strmatch(ref_feature,'phase');
        
        for j = 1:numel(refsig)
            phase(j,:)=mod(angle(refsig{j}),2*pi);
            PHASE(j,:)=phase(j,:);
        end
        
        phase = reshape(phase,numel(refsig),nsamples-(2*round(D.fsample/6))+1, ntrials);  % remove start and end of each trial to avoid filter artefacts
        
        COS(j,:)=((cos(PHASE(j,:))-mean(cos(PHASE(j,:))))./std(cos(PHASE(j,:))))';
        SIN(j,:)=((sin(PHASE(j,:))-mean(sin(PHASE(j,:))))./std(sin(PHASE(j,:))))';
        
        Refsig = cell(1, length(freqoi));
        for j = 1:length(freqoi)
            Refsig{j} = zeros(size(Yr, 1), nsamples-(2*round(D.fsample/6))+1, ntrials);
        end
        
        for i = 1:ntrials
            spectrum = ft_specest_hilbert(squeeze(Yr(:,:, i)), times,...
                'freqoi', freqoi, 'width', lowamp_width, 'filttype', 'but', 'filtorder', 2,...
                'filtdir', 'twopass', 'verbose', 0);
            
            for j = 1:length(freqoi)
                tmp = spm_squeeze(spectrum(:, j, :), 2);
                Refsig{j}(:,:,i) = tmp(:,round(D.fsample/6):end-round(D.fsample/6)); %to remove edge artefacts
            end
        end
        
        for j = 1:numel(Refsig)
            Refsig{j} = reshape(Refsig{j}, size(Refsig{j}, 1), []);
            lowamp(j,:) = abs(Refsig{j});
            LOWAMP(j,:) = ((lowamp(j,:)-mean(lowamp(j,:)))./std(lowamp(j,:)))';
        end
        
        lowamp = reshape(lowamp,numel(refsig),nsamples-(2*round(D.fsample/6))+1, ntrials);
        
    end
    
    J=numel(refsig);
    
end


%%%%%%%%%%%%

switch ref_feature
    case 'amplitude'
        freqoi = S.phasefreq;
        width  = S.phaseres;
        lowamp_width = S.lowampres;
    case 'phase'
        freqoi = S.ampfreq;
        width  = S.ampres;
end

fsample=D.fsample;
fname=D.fname;
clear BF D

%initialize variables
Beta=nan(nphase,namp,nvert,3,ntrials);
r_GLM=nan(nphase,namp,nvert,ntrials);
r_GLM_amp=r_GLM;
r_GLM_total=r_GLM;
all_Beta=nan(nphase,namp,nvert,3);
all_r=nan(nphase,namp,nvert);
all_r_amp=all_r;
all_r_total=all_r;
all_Bnorm=all_r;
pb=nan(nvert,1);
pb_amp=pb;
pb_total=pb;


for i = 1:nvert
    if ~isnan(W{i})
        w    = W{i};
        
        source=w*reshape(Y, nchan, []);
        source=reshape(source,nsamples,ntrials);
        
        try S.reference.within;
            
            if strmatch(ref_feature,'phase')
                
                Foi=S.phasefreq;
                
                for j=1:length(Foi)
                    
                    Yh = 0*source;
                    for i1 = 1:ntrials
                        Yh(: ,i1) = spm_squeeze(ft_specest_hilbert(source(:, i1)', times,...
                            'freqoi', Foi(j), 'width', S.phaseres, 'filttype', 'but', ...
                            'filtorder', 2,  'filtdir', 'twopass', 'verbose', 0), 2);
                        if ismember(i1, Ibar)
                            spm_progress_bar('Set', i1); drawnow;
                        end
                    end
                    
                    phase_tmp = mod(angle(Yh),2*pi);
                    phase(j,:,:) = phase_tmp(round(fsample/6):end-round(fsample/6),:);  % remove start and end of each trial to avoid filter artefacts
                    PHASE(j,:) = reshape(phase(j,:,:),1,[]);
                    
                    Ya = 0*source;
                    for i2 = 1:ntrials
                        Ya(: ,i2) = spm_squeeze(ft_specest_hilbert(source(:, i2)', times,...
                            'freqoi', Foi(j), 'width', S.lowampres, 'filttype', 'but', ...
                            'filtorder', 2,  'filtdir', 'twopass', 'verbose', 0), 2);
                        if ismember(i2, Ibar)
                            spm_progress_bar('Set', i2); drawnow;
                        end
                    end
                    
                    lowamp_tmp = abs(Ya);
                    lowamp(j,:,:) = lowamp_tmp(round(fsample/6):end-round(fsample/6),:);  % remove start and end of each trial to avoid filter artefacts
                    LOWAMP(j,:) = reshape(lowamp(j,:,:),1,[]);
                    LOWAMP(j,:) = ((LOWAMP(j,:)-mean(LOWAMP(j,:)))./std(LOWAMP(j,:)))';
                    COS(j,:)=((cos(PHASE(j,:))-mean(cos(PHASE(j,:))))./std(cos(PHASE(j,:))))';
                    SIN(j,:)=((sin(PHASE(j,:))-mean(sin(PHASE(j,:))))./std(sin(PHASE(j,:))))';
                    
                    J=j;
                end
                
            elseif strmatch(ref_feature,'amplitude')
                
                Foi=S.ampfreq;
                
                for j=1:length(Foi)
                    
                    Yh = 0*source;
                    for i1 = 1:ntrials
                        Yh(: ,i1) = spm_squeeze(ft_specest_hilbert(source(:, i1)', times,...
                            'freqoi', Foi(j), 'width', S.ampres, 'filttype', 'but', ...
                            'filtorder', 2,  'filtdir', 'twopass', 'verbose', 0), 2);
                        if ismember(i1, Ibar)
                            spm_progress_bar('Set', i1); drawnow;
                        end
                    end
                    
                    amplitude_tmp=abs(Yh);
                    amplitude(j,:,:)=amplitude_tmp(round(fsample/6):end-round(fsample/6),:);  % remove start and end of each trial to avoid filter artefacts
                    AMPLITUDE(j,:)=reshape(amplitude(j,:,:),1,[]);
                    AMPLITUDE(j,:)=((AMPLITUDE(j,:)-mean(AMPLITUDE(j,:)))./std(AMPLITUDE(j,:)))';
                    
                    for k=1:ntrials
                        amplitude(j,:,k)=((amplitude(j,:,k)-mean(amplitude(j,:,k)))./std(amplitude(j,:,k)))';
                    end
                end
                J=j;
            end
            
        end
        
        
        for f = 1:length(freqoi)
            
            Yh = 0*source;
            for i1 = 1:ntrials
                Yh(: ,i1) = spm_squeeze(ft_specest_hilbert(source(:, i1)', times,...
                    'freqoi', freqoi(f), 'width', width, 'filttype', 'but', ...
                    'filtorder', 2,  'filtdir', 'twopass', 'verbose', 0), 2);
                if ismember(i1, Ibar)
                    spm_progress_bar('Set', i1); drawnow;
                end
            end
            
            if strmatch(ref_feature,'amplitude')
                
                phase = mod(angle(Yh),2*pi);
                phase = phase(round(fsample/6):end-round(fsample/6),:);  % remove start and end of each trial to avoid filter artefacts
                PHASE = reshape(phase,1,[]);
                
                Ya = 0*source;
                for i2 = 1:ntrials
                    Ya(: ,i2) = spm_squeeze(ft_specest_hilbert(source(:, i2)', times,...
                        'freqoi', freqoi(f), 'width', lowamp_width, 'filttype', 'but', ...
                        'filtorder', 2,  'filtdir', 'twopass', 'verbose', 0), 2);
                    if ismember(i2, Ibar)
                        spm_progress_bar('Set', i2); drawnow;
                    end
                end
                
                lowamp = abs(Ya);
                lowamp = lowamp(round(fsample/6):end-round(fsample/6),:);  % remove start and end of each trial to avoid filter artefacts
                LOWAMP = reshape(lowamp,1,[]);
                LOWAMP = ((LOWAMP-mean(LOWAMP))./std(LOWAMP))';
                COS=((cos(PHASE)-mean(cos(PHASE)))./std(cos(PHASE)))';
                SIN=((sin(PHASE)-mean(sin(PHASE)))./std(sin(PHASE)))';
                
                
            elseif strmatch(ref_feature,'phase')
                
                amplitude=abs(Yh);
                amplitude=amplitude(round(fsample/6):end-round(fsample/6),:);  % remove start and end of each trial to avoid filter artefacts
                AMPLITUDE=reshape(amplitude,1,[]);
                AMPLITUDE=((AMPLITUDE-mean(AMPLITUDE))./std(AMPLITUDE))';
                
                for k=1:ntrials
                    amplitude(:,k)=((amplitude(:,k)-mean(amplitude(:,k)))./std(amplitude(:,k)))';
                end
                
            end
            
            for j = 1:J
                
                for k = 1:ntrials
                    if strmatch(ref_feature,'amplitude')
                        Xk=[(cos(phase(:,k))-mean(cos(phase(:,k)),1))./std(cos(phase(:,k))),(sin(phase(:,k))-mean(sin(phase(:,k)),1))./std(sin(phase(:,k))),(lowamp(:,k)-mean(lowamp(:,k),1))./std(lowamp(:,k))];
                        yk=squeeze(amplitude(j,:,k))';
                    elseif strmatch(ref_feature,'phase')
                        Xk=[(cos(squeeze(phase(j,:,k)))-mean(cos(squeeze(phase(j,:,k)))))./std(cos(squeeze(phase(j,:,k))));(sin(squeeze(phase(j,:,k)))-mean(sin(squeeze(phase(j,:,k)))))./std(sin(squeeze(phase(j,:,k))));(squeeze(lowamp(j,:,k))-mean(squeeze(lowamp(j,:,k))))./std(squeeze(lowamp(j,:,k)))]';
                        yk=squeeze(amplitude(:,k));
                    end
                    
                    %%glm
                    Beta(f,j,i,:,k)=(yk'*pinv(Xk'))';
                    
                end %ntrials
                
                
                %GLM for whole data
                if strmatch(ref_feature,'amplitude')
                    X=[COS,SIN,LOWAMP];
                    y=AMPLITUDE(j,:)';
                elseif strmatch(ref_feature,'phase')
                    X=[COS(j,:);SIN(j,:);LOWAMP(j,:)]';
                    y=AMPLITUDE;
                end
                
                all_Beta(f,j,i,:)=y'*pinv(X');
                
                all_SSy=sum((y-mean(y)).^2);
                all_residuals=y-(all_Beta(f,j,i,1).*X(:,1)+all_Beta(f,j,i,2).*X(:,2));
                all_SSe=sum((all_residuals-mean(all_residuals)).^2);
                all_r(f,j,i)=real(sqrt((all_SSy-all_SSe)/all_SSy));
                
                all_residuals_amp=y-(all_Beta(f,j,i,3).*X(:,3));
                all_SSe_amp=sum((all_residuals_amp-mean(all_residuals_amp)).^2);
                all_r_amp(f,j,i)=real(sqrt((all_SSy-all_SSe_amp)/all_SSy));
                
                all_residuals_total=y-(all_Beta(f,j,i,1).*X(:,1)+all_Beta(f,j,i,2).*X(:,2)+all_Beta(f,j,i,3).*X(:,3));
                all_SSe_total=sum((all_residuals_total-mean(all_residuals_total)).^2);
                all_r_total(f,j,i)=real(sqrt((all_SSy-all_SSe_total)/all_SSy));
                
                all_Bnorm(f,j,i)=sqrt(all_Beta(f,j,i,1).^2+all_Beta(f,j,i,2).^2);
                
            end %numel(refsig)
        end %~isnan(W{i})
        
    end %nvert
end %freqoi

% avg over spectrum
IM_all_r=squeeze(mean(mean(all_r,2),1));
IM_all_r_amp=squeeze(mean(mean(all_r_amp,2),1));
IM_all_r_total=squeeze(mean(mean(all_r_total,2),1));
IM_all_B3=squeeze(mean(mean(all_Beta(:,:,:,3),2),1));
IM_all_B2=squeeze(mean(mean(all_Beta(:,:,:,2),2),1));
IM_all_B1=squeeze(mean(mean(all_Beta(:,:,:,1),2),1));

IM_trials_B3=squeeze(mean(mean(Beta(:,:,:,3,:),2),1));
IM_trials_B2=squeeze(mean(mean(Beta(:,:,:,2,:),2),1));
IM_trials_B1=squeeze(mean(mean(Beta(:,:,:,1,:),2),1));

B1=squeeze(mean(mean(Beta(:,:,:,1,:),2),1));
B2=squeeze(mean(mean(Beta(:,:,:,2,:),2),1));
B3=squeeze(mean(mean(Beta(:,:,:,3,:),2),1));

% seond level stats
for i=1:nvert
    
    V=[];
    Xb=[];
    Xb(1:ntrials,1)=ones(ntrials,1);Xb(ntrials+1:2*ntrials,2)=ones(ntrials,1);
    yb=[B1(i,:),B2(i,:)];
    c=eye(2);
    [F,df,Beta_b,xX,xCon]=spm_ancova(Xb,V,yb',c);
    pb(i)=1-spm_Fcdf(F,df(1),df(2));
    
    %to test for significance amp-amp
    yb=[B3(i,:)];
    [H,P] = ttest(yb);
    pb_amp(i)=P;
    
    %to test for significance phase-amp & amp-amp
    Xb=[];
    Xb(1:ntrials,1)=ones(ntrials,1);Xb(ntrials+1:2*ntrials,2)=ones(ntrials,1);Xb(2*ntrials+1:3*ntrials,3)=ones(ntrials,1);
    c=eye(3);
    yb=[B1(i,:),B2(i,:),B3(i,:)];
    [F_total,df,Beta_b,xX,xCon]=spm_ancova(Xb,V,yb',c);
    pb_total(i)=1-spm_Fcdf(F_total,df(1),df(2));
    
end


%%write out images

cnt=1;

image(cnt).val     = pb;
image(cnt).label   = ['p_pac_'  spm_file(fname, 'basename')];
cnt=cnt+1;
image(cnt).val     = pb_amp;
image(cnt).label   = ['p_amp_'  spm_file(fname, 'basename')];
cnt=cnt+1;
image(cnt).val     = pb_total;
image(cnt).label   = ['p_total_'  spm_file(fname, 'basename')];
cnt=cnt+1;
image(cnt).val     = IM_all_r;
image(cnt).label   = ['r_pac_'  spm_file(fname, 'basename')];
cnt=cnt+1;
image(cnt).val     = IM_all_r_amp;
image(cnt).label   = ['r_amp_'  spm_file(fname, 'basename')];
cnt=cnt+1;
image(cnt).val     = IM_all_r_total;
image(cnt).label   = ['r_total_'  spm_file(fname, 'basename')];
cnt=cnt+1;
image(cnt).val     = IM_all_B3;
image(cnt).label   = ['B3_'  spm_file(fname, 'basename')];
cnt=cnt+1;
image(cnt).val     = IM_all_B2;
image(cnt).label   = ['B2_'  spm_file(fname, 'basename')];
cnt=cnt+1;
image(cnt).val      = IM_all_B1;
image(cnt).label   = ['B1_'  spm_file(fname, 'basename')];
cnt=cnt+1;

% for k=1:ntrials
%     image(cnt).val = IM_trials_B3(:,k);
%     image(cnt).label   = ['trial',num2str(k),'_B3_'  spm_file(fname, 'basename')];
%     cnt=cnt+1;
%     image(cnt).val = IM_trials_B2(:,k);
%     image(cnt).label   = ['trial',num2str(k),'_B2_'  spm_file(fname, 'basename')];
%     cnt=cnt+1;
%     image(cnt).val = IM_trials_B1(:,k);
%     image(cnt).label   = ['trial',num2str(k),'_B1_'  spm_file(fname, 'basename')];
%     cnt=cnt+1;
% end


res = image;