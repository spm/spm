function res = bf_output_image_pac(BF, S)
% Computes phase-amplitude coupling
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Bernadette van Wijk, Vladimir Litvak
% $Id: bf_output_image_pac.m 7703 2019-11-22 12:06:29Z guillaume $

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
    
    shuffle         = cfg_entry;
    shuffle.tag     = 'shuffle';
    shuffle.name    = 'Shuffle';
    shuffle.strtype = 'w';
    shuffle.num     = [1 1];
    shuffle.help    = {'Shuffle the reference channel to produce the null case.',...
        'Specify the number of shufflings'};
    shuffle.val = {0};
    
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
    
    image_pac      = cfg_branch;
    image_pac.tag  = 'image_pac';
    image_pac.name = 'PAC image';
    image_pac.val  = {whatconditions, sametrials, shuffle, woi, phasefreq, ....
        phaseres, ampfreq, ampres, reference, modality, ncomp_amp, ncomp_phase, outputname}; 
    
    res = image_pac;
    
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
    if nphase>=namp
        ref_feature = 'amplitude';
    else
        ref_feature = 'phase';
    end
    
    Yr = Y;
end

switch ref_feature
    case 'amplitude'
        freqoi = S.ampfreq;
        width  = S.ampres;
    case 'phase'
        freqoi = S.phasefreq;
        width  = S.phaseres;
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

switch ref_feature
    case 'amplitude'
        freqoi = S.phasefreq;
        width  = S.phaseres;
    case 'phase'
        freqoi = S.ampfreq;
        width  = S.ampres;
end

if ncomponents>1
    pac = nan(nphase, namp, nvert, ncomponents);
else
    pac = nan(nphase, namp, nvert);
end

if S.shuffle
    if ncomponents>1
        spac = repmat(pac, [1 1 1 1 S.shuffle]);
    else
        spac = repmat(pac, [1 1 1 S.shuffle]);
    end
    sind = zeros(S.shuffle, ntrials);
    for s = 1:S.shuffle
        sind(s, :) = randperm(ntrials);
    end
end

for f = 1:length(freqoi)
    fprintf('%s of %s\n',num2str(freqoi(f)),num2str(freqoi(end)));
    
    spm_progress_bar('Init', ntrials, ...
        sprintf('Computing data spectra')); drawnow;
    if ntrials > 100, Ibar = floor(linspace(1, ntrials,100));
    else Ibar = 1:ntrials; end
    
    
    Yh = 0*Y;
    for i = 1:ntrials
        Yh(: , : ,i) = spm_squeeze(ft_specest_hilbert(squeeze(Y(:,:, i)), times,...
            'freqoi', freqoi(f), 'width', width, 'filttype', 'but', ...
            'filtorder', 2,  'filtdir', 'twopass', 'verbose', 0), 2);
        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end
    end
    
    Yh = Yh(:,round(D.fsample/6):end-round(D.fsample/6),:);  % remove start and end of each trial to avoid filter artefacts
    Yh = reshape(Yh, nchan, []);
    
    spm_progress_bar('Clear');
    
    spm_progress_bar('Init', nvert, ...
        sprintf('Scanning grid points image')); drawnow;
    if nvert > 100, Ibar = floor(linspace(1, nvert,100));
    else Ibar = 1:nvert; end
    
    for i = 1:nvert
        if ~isnan(W{i})
            w    = W{i};
            
            sYh  = w*Yh;
            
            for j = 1:numel(refsig)
                if size(refsig{j}, 1)>1
                    rYh = w*refsig{j};
                else
                    rYh = refsig{j};
                end
                
                for shuffle = 0:S.shuffle
                    if shuffle
                        nc=size(rYh,1);
                        rYh = reshape(rYh, nc, nsamples-(2*round(D.fsample/6))+1, ntrials);
                        rYh = rYh(:,:, sind(shuffle, :));
                        rYh = reshape(rYh,nc,[]);
                    end
                    
                    switch ref_feature
                        case 'amplitude'
                            amp = abs(rYh);
                            phase = mod(angle(sYh),2*pi);
                        case 'phase'
                            amp = abs(sYh);
                            phase = mod(angle(rYh),2*pi);
                    end
                    
                    if S.ncomp_amp==1 && S.ncomp_phase==1;
                        amp(1,:)=amp(1,:)-mean(amp(1,:)); %(1,:) so it also works when keeping more components but computing PAC for just 1
                        cpac = abs(sum(amp(1,:).*exp(sqrt(-1)*phase(1,:)))/(nsamples*ntrials));
                        cpac = cpac*sqrt(2/var(amp(1,:)));
                    end
                    if S.ncomp_amp==1 && S.ncomp_phase>1;
                        amp(1,:)=amp(1,:)-mean(amp(1,:));
                        for kk=1:S.ncomp_phase
                            cpac(kk) = abs(sum(amp(1,:).*exp(sqrt(-1)*phase(kk,:)))/(nsamples*ntrials));
                            cpac(kk) = cpac(kk)*sqrt(2/var(amp(1,:)));
                        end
                    end
                    if S.ncomp_amp>1 && S.ncomp_phase==1;
                        for kk=1:S.ncomp_amp
                            amp(kk,:)=amp(kk,:)-mean(amp(kk,:));
                            cpac(kk) = abs(sum(amp(kk,:).*exp(sqrt(-1)*phase(1,:)))/(nsamples*ntrials));
                            cpac(kk) = cpac(kk)*sqrt(2/var(amp(kk,:)));
                        end
                    end
                    if S.ncomp_amp>1 && S.ncomp_phase>1;
                        cnt=1;
                        for k1=1:S.ncomp_amp
                            for k2=1:S.ncomp_phase
                                cpac(cnt) = abs(sum(amp(k1,:).*exp(sqrt(-1)*phase(k2,:)))/(nsamples*ntrials));
                                cpac(cnt) = cpac(cnt)*sqrt(2/var(amp(k1,:)));
                                cnt=cnt+1;
                            end
                        end
                    end
                    
                    
                    switch ref_feature
                        case 'amplitude'
                            if shuffle
                                if ncomponents>1;
                                    spac(f, j, i, :, shuffle) = cpac;
                                else
                                    spac(f, j, i, shuffle) = cpac;
                                end
                            else
                                if ncomponents>1;
                                    pac(f,j,i,:)=cpac;
                                else
                                    pac(f, j, i) = cpac;
                                end
                            end
                        case 'phase'
                            if shuffle
                                if ncomponents>1;
                                    spac(j, f, i, :, shuffle) = cpac;
                                else
                                    spac(j, f, i, shuffle) = cpac;
                                end
                            else
                                if ncomponents>1;
                                    pac(j, f, i, :) = cpac;
                                else
                                    pac(j, f, i) = cpac;
                                end
                            end
                    end
                end
            end
        end
        
        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end
    end
    
    spm_progress_bar('Clear');
end


c = 1;
for kk=1:ncomponents
    
    if max(nphase, namp)>1
        if ncomponents>1
            image(c).val     = squeeze(sum(sum(pac(:,:,:,kk), 2), 1));
            image(c).label   = ['pac_total_MEGcomponent' num2str(kk) '_' S.outputname spm_file(D.fname, 'basename')];
        else
            image(c).val     = squeeze(sum(sum(pac, 2), 1));
            image(c).label   = ['pac_total_'  spm_file(D.fname, 'basename')];
        end
        c = c+1;
    end
    
    for f = 1:nphase
        for g = 1:namp
            if ncomponents>1
                image(c).val     = squeeze(pac(f, g, :,kk));
                image(c).label   = ['pac_phase_' num2str(S.phasefreq(f)) 'Hz_amp_' num2str(S.ampfreq(g)) 'Hz_MEGcomponent' num2str(kk) '_' S.outputname spm_file(D.fname, 'basename')];
            else
                image(c).val     = squeeze(pac(f, g, :));
                image(c).label   = ['pac_phase_' num2str(S.phasefreq(f)) 'Hz_amp_' num2str(S.ampfreq(g)) 'Hz_' S.outputname spm_file(D.fname, 'basename')];
            end
            c = c+1;
        end
    end
    
    for shuffle = 1:S.shuffle
        
        if max(nphase, namp)>1
            if ncomponents>1
                image(c).val     = squeeze(sum(sum(spac(:,:,:,kk,shuffle), 2), 1));
                image(c).label   = ['shuffled' num2str(shuffle) '_pac_total_MEGcomponent' num2str(kk) '_' spm_file(D.fname, 'basename')];
            else
                image(c).val     = squeeze(sum(sum(spac(:,:,:,shuffle), 2), 1));
                image(c).label   = ['shuffled' num2str(shuffle) '_pac_total_'  spm_file(D.fname, 'basename')];
            end
            c = c+1;
        end
        
        for f = 1:nphase
            for g = 1:namp
                if ncomponents>1
                    image(c).val     = squeeze(spac(f, g, :, kk, shuffle));
                    image(c).label   = ['shuffled' num2str(shuffle) '_pac_phase_' num2str(S.phasefreq(f)) 'Hz_amp_' num2str(S.ampfreq(g)) 'Hz_MEGcomponent' num2str(kk) '_' S.outputname spm_file(D.fname, 'basename')];
                else
                    image(c).val     = squeeze(spac(f, g, :, shuffle));
                    image(c).label   = ['shuffled' num2str(shuffle) '_pac_phase_' num2str(S.phasefreq(f)) 'Hz_amp_' num2str(S.ampfreq(g)) 'Hz_' S.outputname spm_file(D.fname, 'basename')];
                end
                c = c+1;
            end
        end
    end
end

res = image;
