function res = bf_write_spmeeg(BF, S)
% Writes out beamformer results as M/EEG dataset
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_write_spmeeg.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0
    mode         = cfg_menu;
    mode.tag     = 'mode';
    mode.name    = 'Writing mode';
    mode.help    = {'How to generate the output'};
    mode.labels  = {
        'Write new'
        'Online montage original data'
        'Copy + online montage'
        }';
    mode.values  = {
        'write'
        'online'
        'onlinecopy'
        }';
    mode.val = {'write'};
    
    modality         = cfg_menu;
    modality.tag     = 'modality';
    modality.name    = 'Modality';
    modality.help    = {'What modality to output'};
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
    
    none = cfg_const;
    none.tag = 'none';
    none.name = 'None';
    none.val  = {0};
    
    addchannels      = cfg_choice;
    addchannels.tag  = 'addchannels';
    addchannels.name = 'Extra channels to add';
    addchannels.values  = {none, spm_cfg_eeg_channel_selector};
    addchannels.val  = {none};
    
    prefix         = cfg_entry;
    prefix.tag     = 'prefix';
    prefix.name    = 'Filename Prefix';
    prefix.help    = {'Specify the string to be prepended to the output (if relevant).'};
    prefix.strtype = 's';
    prefix.num     = [1 Inf];
    prefix.val     = {'B'};
    
    spmeeg      = cfg_branch;
    spmeeg.tag  = 'spmeeg';
    spmeeg.name = 'SPM M/EEG dataset';
    spmeeg.val  = {mode, modality, addchannels, prefix};
    
    res = spmeeg;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;

if isfield(S.addchannels, 'channels')
    addchannels = D.chanlabels(D.selectchannels(spm_cfg_eeg_channel_selector(S.addchannels.channels)));
else
    addchannels = {};
end

usemontage = 0; 

if isfield(BF.output, 'montage')
    
    usemontage = 1;
    
    montage = BF.output.montage.(S.modality);
    
    [sel1, sel2] = spm_match_str(montage.labelorg, addchannels);
    for i = 1:length(sel1)
        montage.labelnew(end+1) = addchannels(sel2(i));
        montage.tra(end+1, sel1(i)) = 1;
        
        montage.chantypenew(end+1) = chantype(D, D.indchannel(addchannels(sel2(i))))';
        montage.chanunitnew(end+1) = units(D, D.indchannel(addchannels(sel2(i))))';

        
    end
    addchannels(sel2) = [];
    
    if ~isempty(addchannels)
        montage.labelorg = [montage.labelorg(:); addchannels(:)];
        montage.labelnew = [montage.labelnew(:); addchannels(:)];
        
        montage.chantypenew = [montage.chantypenew(:); chantype(D, D.indchannel(addchannels))'];
        montage.chanunitnew = [montage.chanunitnew(:); units(D, D.indchannel(addchannels))'];
        
        na = numel(addchannels);
        montage.tra((end+1):(end+na), (end+1):(end+na)) = eye(na);
    end
    
    montage.chantypeold = chantype(D, D.indchannel(montage.labelorg))';
    montage.chanunitold = units(D, D.indchannel(montage.labelorg))';
elseif isfield(BF.output, 'sourcedata')
    ftdata = BF.output.sourcedata.(S.modality).ftdata;
    
    if ~isempty(addchannels)
        addind = D.indchannel(addchannels);
        ftdata.label = [ftdata.label;D.chanlabels(addind)'];
        if numel(ftdata.trial) == length(BF.features.trials)
            for i = 1:numel(ftdata.trial)
                ftdata.trial{i} = [ftdata.trial{i};D(addind, D.indsample(ftdata.time{i}(1)):D.indsample(ftdata.time{i}(end)), BF.features.trials(i))];
            end
        else
            error('Cannot match trials for adding extra channels.')
        end
    end
    
    Ds = spm_eeg_ft2spm(ftdata, [S.prefix D.fname]);
    Ds = trialonset(Ds, ':', D.trialonset);
    
    if isfield(BF.output.sourcedata.(S.modality), 'events')
        for i = 1:Ds.ntrials
            evold = events(D, BF.features.trials(i));
            if iscell(evold)
                evold = evold{1};
            end
            
            evnew = BF.output.sourcedata.(S.modality).events;
            if iscell(evnew)
                evnew = evnew{i};
            end
            
            ev = spm_cat_struct(evold, evnew);
            
            Ds = events(Ds, i, ev);
        end
    else
        if ~isempty(D.events)
            Ds = events(Ds, ':', D.events);
        end
    end
    
    
    if isfield(BF.sources, 'voi')
        Ds = chantype(Ds, ':', 'LFP');
    else
        Ds = chantype(Ds, ':', 'SRC');
    end
    
    if ~isempty(addchannels)
         Ds = chantype(Ds, (Ds.nchannels-length(addind)+1):Ds.nchannels, D.chantype(addind));
    end
    
    save(Ds);
    
    D = Ds;
else
    error('Unsupported option');
end

if usemontage
    if D.montage('getindex')
        vmontage = rmfield(D.montage('getmontage'), 'channels');
        montage  =  ft_apply_montage(vmontage, montage);
    end
    
    S1 = [];
    S1.montage = montage;
    S1.prefix  = S.prefix; % ignored for online
    S1.keepsensors = false;
    S1.keepothers  = false;
        
    switch S.mode
        case 'write'                                              
            S1.mode    = 'write';
        case 'online'
            S1.mode = 'switch';
        case 'onlinecopy' 
            S1.mode    = 'switch';            
            D = copy(D, [S.prefix D.fname]);
    end        
    
    S1.D = D;
    D = spm_eeg_montage(S1);
end

res.files = {fullfile(D)};
