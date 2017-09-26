function D = spm_eeg_prep(S)
% Prepare converted M/EEG data for further analysis
% FORMAT D = spm_eeg_prep(S)
% S                 - configuration structure (optional)
% (optional) fields of S:
%   S.D             - MEEG object or filename of M/EEG mat-file
%   S.task          - action string. One of 'settype', 'defaulttype',
%                     'loadtemplate','setcoor2d', 'project3d', 'loadeegsens',
%                     'defaulteegsens', 'sens2chan', 'headshape',
%                     'coregister', 'sortconditions'
%
%   S.updatehistory - update history information [default: true]
%   S.save          - save MEEG object [default: false]
%
% D                 - MEEG object
%__________________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: spm_eeg_prep.m 7175 2017-09-26 14:50:54Z vladimir $

D = spm_eeg_load(S.D);

switch lower(S.task)
    %----------------------------------------------------------------------
    case 'setbadchan'
        %----------------------------------------------------------------------
        D = badchannels(D, D.selectchannels(S.channels), S.status);
        %----------------------------------------------------------------------
    case 'settype'
        %----------------------------------------------------------------------
        D = chantype(D, S.ind, S.type);
        
        %----------------------------------------------------------------------
    case 'defaulttype'
        %----------------------------------------------------------------------
        if isfield(S, 'ind')
            ind = S.ind;
        else
            ind = 1:D.nchannels;
        end
        
        dictionary = {
            'eog',           'EOG';
            'eeg',           'EEG';
            'ecg',           'ECG';
            'lfp',           'LFP';
            'emg',           'EMG';
            'meg',           'MEG';
            'ref',           'REF';
            'megref'         'REF';
            'megmag',        'MEGMAG';
            'megplanar',     'MEGPLANAR';
            'meggrad',       'MEGGRAD';
            'refmag',        'REFMAG';
            'refgrad',       'REFGRAD'
            'refplanar'      'REFPLANAR'
            };
        
        D = chantype(D, ind, 'Other');
        
        type = ft_chantype(D.chanlabels);
        
        % If there is useful information in the original types it
        % overwrites the default assignment
        if isfield(D, 'origchantypes')
            [sel1, sel2] = spm_match_str(chanlabels(D, ind), D.origchantypes.label);
            
            type(ind(sel1)) = D.origchantypes.type(sel2);
        end
        
        spmtype = repmat({'Other'}, 1, length(ind));
        
        [sel1, sel2] = spm_match_str(type(ind), dictionary(:, 1));
        
        spmtype(sel1) = dictionary(sel2, 2);
        
        D = chantype(D, ind, spmtype);
        
        %----------------------------------------------------------------------
        
    case 'bidschantype'
        
        dictionary = { % Should be updates with updates to BIDS specification
            'MEGMAG'   'MEGMAG'
            'MEGGRAD'  'MEGPLANAR'
            'MEG'      'MEG'
            'MEGREF'   'REF'
            'EEG'      'EEG'
            'EOG'      'EOG'
            'ECG'      'ECG'
            'EMG'      'EMG'
            'TRIG'     'Other'
            'AUDIO'    'Other'
            'PD'       'Other'
            'ET'       'Other'
            'MISC'     'Other'
            };
        
        bids_chan = spm_load(S.filename);
        
        [sel1, sel2] = spm_match_str(D.chanlabels, bids_chan.name);
        
        type = bids_chan.type(sel2);
                              
        [sel3, sel4] = spm_match_str(type, dictionary(:, 1));
        
        type(sel3) = dictionary(sel4, 2);
        
        D = chantype(D, sel1, type);  
        
    case 'bidschanstatus'
        bids_chan = spm_load(S.filename);
        
        [sel1, sel2] = spm_match_str(D.chanlabels, bids_chan.name);
        
        sel3 = strmatch('bad', bids_chan.status(sel2));
        
        D = badchannels(D, sel1, 0);
        
        D = badchannels(D, sel1(sel3), 1);
        
    case {'loadtemplate', 'setcoor2d', 'project3d'}
        %----------------------------------------------------------------------
        chanind = 1:D.nchannels;
        
        switch lower(S.task)
            case 'loadtemplate'
                template    = load(S.P); % must contain Cpos, Cnames
                xy          = template.Cpos;
                label       = template.Cnames;
            case 'setcoor2d'
                xy          = S.xy;
                label       = S.label;
            case 'project3d'
                if ~isfield(D, 'val')
                    D.val = 1;
                end
                if isfield(D, 'inv') && isfield(D.inv{D.val}, 'datareg')
                    datareg = D.inv{D.val}.datareg;
                    ind     = strmatch(S.modality, {datareg(:).modality}, 'exact');
                    sens    = datareg(ind).sensors;
                else
                    sens    = D.sensors(S.modality);
                end
                [xy, label] = spm_eeg_project3D(sens, S.modality);
                
                switch S.modality
                    case 'EEG'
                        chanind  = D.indchantype('EEG');
                    case 'MEG'
                        chanind  = D.indchantype('MEGANY');
                    case 'MEGCOMB'
                        chanind  = D.indchantype('MEGCOMB');
                end
        end
        
        megcombind   = D.indchantype('MEGCOMB');
        
        if ~isempty(megcombind)
            chanset = spm_eeg_planarchannelset(D.chanlabels);
            [sel1, sel2] = spm_match_str(lower(chanset(:, 1)), lower(label));
            [sel3, sel4] = spm_match_str(lower(chanset(:, 2)), lower(label));
            
            [sel5, sel6, sel7] = intersect(sel1, sel3);
            
            label = [label(:); chanset(sel5, 3)];
            
            xy    = [xy 0.5*(xy(:, sel2(sel6)) + xy(:, sel4(sel7)))];
        end
        
        [sel1, sel2] = spm_match_str(lower(D.chanlabels(chanind)), lower(label));
        sel1         = chanind(sel1);        
        
        if ~isempty(sel1)
            
            megind = D.indchantype('MEG');
            eegind = D.indchantype('EEG');
            
            if ~isempty(intersect(megind, sel1)) && ~isempty(setdiff(megind, sel1))
                error('2D locations not found for all MEG channels');
            end
            
            if ~isempty(intersect(megcombind, sel1)) && ~isempty(setdiff(megcombind, sel1))
                error('2D locations not found for all MEGCOMB channels');
            end
            
            if ~isempty(intersect(eegind, sel1)) && ~isempty(setdiff(eegind, sel1))
                warning(['2D locations not found for all EEG channels, changing type of channels', ...
                    num2str(setdiff(eegind, sel1)) ' to ''Other''']);
                
                D = chantype(D, setdiff(eegind, sel1), 'Other');
            end
            
            D = coor2D(D, sel1, num2cell(xy(:, sel2)));
        end
        
        %----------------------------------------------------------------------
    case 'loadeegsens'
        %----------------------------------------------------------------------
        switch S.source
            case 'mat'
                senspos = load(S.sensfile);
                name    = fieldnames(senspos);
                senspos = getfield(senspos,name{1});
                
                label = chanlabels(D, D.indchantype('EEG'));
                
                if size(senspos, 1) ~= length(label)
                    error('To read sensor positions without labels the numbers of sensors and EEG channels should match.');
                end
                
                elec = [];
                elec.chanpos = senspos;
                elec.elecpos = senspos;
                elec.label = label;
                
                headshape = load(S.headshapefile);
                name    = fieldnames(headshape);
                headshape = getfield(headshape,name{1});
                
                shape = [];
                
                fidnum = 0;
                while ~all(isspace(S.fidlabel))
                    fidnum = fidnum+1;
                    [shape.fid.label{fidnum},S.fidlabel] = strtok(S.fidlabel);
                end
                
                if (fidnum < 3)  || (size(headshape, 1) < fidnum)
                    error('At least 3 labeled fiducials are necessary');
                end
                
                shape.fid.pnt = headshape(1:fidnum, :);
                
                if size(headshape, 1) > fidnum
                    shape.pnt = headshape((fidnum+1):end, :);
                else
                    shape.pnt = [];
                end
            case 'locfile'
                label = chanlabels(D, D.indchantype('EEG'));
                
                elec = ft_read_sens(S.sensfile);
                
                % Remove headshape points
                hspind = strmatch('headshape', elec.label);
                elec.chanpos(hspind, :) = [];
                elec.elecpos(hspind, :) = [];
                elec.chantype(hspind, :) = [];
                elec.chanunit(hspind, :) = [];
                elec.label(hspind)  = [];
                
                % This handles FIL Polhemus case and other possible cases
                % when no proper labels are available.
                if isempty(intersect(label, elec.label))
                    ind = str2num(strvcat(elec.label));
                    if length(ind) == length(label)
                        elec.label = label(ind);
                    else
                        error('To read sensor positions without labels the numbers of sensors and EEG channels should match.');
                    end
                end
                
                shape = spm_eeg_fixpnt(ft_read_headshape(S.sensfile));
                
                % In case electrode file is used for fiducials, the
                % electrodes can be used as headshape
                if ~isfield(shape, 'pnt') || isempty(shape.pnt) && ...
                        size(shape.fid.pnt, 1) > 3
                    shape.pnt = shape.fid.pnt;
                end
                
        end
        
        elec = ft_convert_units(elec, 'mm');
        shape= ft_convert_units(shape, 'mm');
        
        if isequal(D.modality(1, 0), 'Multimodal')
            if ~isempty(D.fiducials) && isfield(S, 'regfid') && ~isempty(S.regfid)
                M1 = coreg(D.fiducials, shape, S.regfid);
                elec = ft_transform_sens(M1, elec);
            else
                error(['MEG fiducials matched to EEG fiducials are required '...
                    'to add EEG sensors to a multimodal dataset.']);
            end
        else
            D = fiducials(D, shape);
        end
        
        D = sensors(D, 'EEG', elec);
        
        %----------------------------------------------------------------------
    case 'defaulteegsens'
        %----------------------------------------------------------------------
        
        template_sfp = dir(fullfile(spm('dir'), 'EEGtemplates', '*.sfp'));
        template_sfp = {template_sfp.name};
        
        if ft_senstype(D.chanlabels(D.indchantype('EEG')), 'ext1020')
            ind = strmatch('ext1020.sfp', template_sfp, 'exact');
        else
            ind = strmatch([ft_senstype(D.chanlabels(D.indchantype('EEG'))) '.sfp'], template_sfp, 'exact');
        end
        
        if ~isempty(ind)
            fid = D.fiducials;
            
            if isequal(D.modality(1, 0), 'Multimodal') && ~isempty(fid)
                
                nzlbl = {'fidnz', 'nz', 'nas', 'nasion', 'spmnas'};
                lelbl = {'fidle', 'fidt9', 'lpa', 'lear', 'earl', 'le', 'l', 't9', 'spmlpa'};
                relbl = {'fidre', 'fidt10', 'rpa', 'rear', 'earr', 're', 'r', 't10', 'spmrpa'};
                
                [sel1, nzind] = spm_match_str(nzlbl, lower(fid.fid.label));
                if ~isempty(nzind)
                    nzind = nzind(1);
                end
                [sel1, leind] = spm_match_str(lelbl, lower(fid.fid.label));
                if ~isempty(leind)
                    leind = leind(1);
                end
                [sel1, reind] = spm_match_str(relbl, lower(fid.fid.label));
                if ~isempty(reind)
                    reind = reind(1);
                end
                
                regfid = fid.fid.label([nzind, leind, reind]);
                if numel(regfid) < 3
                    error('Could not automatically understand the MEG fiducial labels. Please use the GUI.');
                else
                    regfid = [regfid(:) {'spmnas'; 'spmlpa'; 'spmrpa'}];
                end
                
                S1 = [];
                S1.D = D;
                S1.task = 'loadeegsens';
                S1.source = 'locfile';
                S1.regfid = regfid;
                S1.sensfile = fullfile(spm('dir'), 'EEGtemplates', template_sfp{ind});
                S1.updatehistory = 0;
                D = spm_eeg_prep(S1);
            else
                elec = ft_read_sens(fullfile(spm('dir'), 'EEGtemplates', template_sfp{ind}));
                
                [sel1, sel2] = spm_match_str(lower(D.chanlabels), lower(elec.label));
                
                sens = elec;
                sens.chanpos  = sens.chanpos(sel2, :);
                sens.elecpos  = sens.elecpos(sel2, :);
                sens.chantype = sens.chantype(sel2, :);
                sens.chanunit = sens.chanunit(sel2, :);
                
                % This takes care of possible case mismatch
                sens.label = D.chanlabels(sel1);
                
                sens.label = sens.label(:);
                
                D = sensors(D, 'EEG', sens);
                
                % Assumes that the first 3 points in standard location files
                % are the 3 fiducials (nas, lpa, rpa)
                fid = [];
                fid.pnt = elec.elecpos;
                fid.fid.pnt = elec.elecpos(1:3, :);
                fid.fid.label = elec.label(1:3);
                
                [xy, label] = spm_eeg_project3D(D.sensors('EEG'), 'EEG');
                
                [sel1, sel2] = spm_match_str(lower(D.chanlabels), lower(label));
                
                if ~isempty(sel1)
                    
                    eegind = strmatch('EEG', chantype(D), 'exact');
                    
                    if ~isempty(intersect(eegind, sel1)) && ~isempty(setdiff(eegind, sel1))
                        warning(['2D locations not found for all EEG channels, changing type of channels ', ...
                            num2str(setdiff(eegind(:)', sel1(:)')) ' to ''Other''']);
                        
                        D = chantype(D, setdiff(eegind, sel1), 'Other');
                    end
                    
                    if any(any(coor2D(D, sel1) - xy(:, sel2)))
                        D = coor2D(D, sel1, num2cell(xy(:, sel2)));
                    end
                end
                
                if strcmp(D.modality(1, 0), 'Multimodal') && ~isempty(D.fiducials)...
                        && isfield(S, 'regfid') && ~isempty(S.regfid)
                    M1 = coreg(D.fiducials, fid, S.regfid);
                    D = sensors(D, 'EEG', ft_transform_sens(M1, D.sensors('EEG')));
                else
                    D = fiducials(D, fid);
                end
            end
        end
        
        %----------------------------------------------------------------------
    case 'loadmegsens'
        %----------------------------------------------------------------------
        hdr  = ft_read_header(S.source);
        D = sensors(D, 'MEG', hdr.grad);
        D = fiducials(D, ft_convert_units(ft_read_headshape(S.source), 'mm'));
        
        if ~isempty(D.indchantype('MEG')) && ~isempty(D.sensors('MEG'))
            
            S1 = [];
            S1.task = 'project3D';
            S1.modality = 'MEG';
            S1.updatehistory = 0;
            S1.D = D;
            
            D = spm_eeg_prep(S1);
        end
        
        %----------------------------------------------------------------------
    case 'sens2chan'
        %----------------------------------------------------------------------
        if isfield(S, 'montage')
            montage = S.montage;
            if ischar(montage)
                montage = getfield(load(montage), 'montage');
            end
        elseif isfield(S, 'refelec');
            sens = sensors(D, 'EEG');
            if isempty(sens)
                error('The montage cannod be applied - no EEG sensors specified');
            end
            if ismember('all', S.refelec)
                refind = 1:numel(sens.label);
            else
                refind = spm_match_str(sens.label, S.refelec);
            end
            
            tra            = eye(numel(sens.label));
            tra(:, refind) = tra(:, refind) - 1/length(refind);
            
            montage.tra = tra;
            montage.labelorg = sens.label;
            montage.labelnew = sens.label;
        else
            error('Montage or list of reference sensors should be specified');
        end
        
        modalities = {'EEG', 'MEG'};
        for m = 1:numel(modalities)
            sens = sensors(D, modalities{m});
            if ~isempty(sens) && ~isempty(intersect(sens.label, montage.labelorg))
                sens = sensors(D, 'EEG');
                sens = ft_apply_montage(sens, montage, 'keepunused', 'no');
                D = sensors(D, modalities{m}, sens);
            end
        end
        
        %----------------------------------------------------------------------
    case 'headshape'
        %----------------------------------------------------------------------
        switch S.source
            case 'mat'
                headshape = load(S.headshapefile);
                name    = fieldnames(headshape);
                headshape = getfield(headshape,name{1});
                
                shape = [];
                
                fidnum = 0;
                while ~all(isspace(S.fidlabel))
                    fidnum = fidnum+1;
                    [shape.fid.label{fidnum},S.fidlabel] = strtok(S.fidlabel);
                end
                
                if (fidnum < 3)  || (size(headshape, 1) < fidnum)
                    error('At least 3 labeled fiducials are necessary');
                end
                
                shape.fid.pnt = headshape(1:fidnum, :);
                
                if size(headshape, 1) > fidnum
                    shape.pnt = headshape((fidnum+1):end, :);
                else
                    shape.pnt = [];
                end
            otherwise
                shape = spm_eeg_fixpnt(ft_read_headshape(S.headshapefile));
                
                % In case electrode file is used for fiducials, the
                % electrodes can be used as headshape
                if ~isfield(shape, 'pnt') || isempty(shape.pnt) && ...
                        size(shape.fid.pnt, 1) > 3
                    shape.pnt = shape.fid.pnt;
                end
        end
        
        shape = ft_convert_units(shape, 'mm');
        
        fid = D.fiducials;
        
        if ~isempty(fid) && isfield(S, 'regfid') && ~isempty(S.regfid)
            M1 = coreg(fid, shape, S.regfid);
            shape = ft_transform_headshape(M1, shape);
        end
        
        D = fiducials(D, shape);
        %----------------------------------------------------------------------
    case 'coregister'
        %----------------------------------------------------------------------
        [D, ok] = check(D, '3d');
        
        if ~ok
            error('Coregistration cannot be performed due to missing data');
        end
        
        try
            val = D.val;
            Msize = D.inv{val}.mesh.Msize;
        catch
            val = 1;
            Msize = 1;
        end
        
        D = spm_eeg_inv_mesh_ui(D, val, 1, Msize);
        D = spm_eeg_inv_datareg_ui(D, val);
        
    case 'sortconditions'
        if ischar(S.condlist)
            cl = getfield(load(S.condlist), 'condlist');
        else
            cl = S.condlist;
        end
        
        D = condlist(D, cl);
        %----------------------------------------------------------------------
    case 'loadbidsevents'
        if ~isequal(D.type, 'continuous')
            error('This operation can only be applied to continuous datasets');
        end
        
        ev_bids = spm_load(S.filename);
        
        ev_spm = struct('type', repmat({'BIDS'}, length(ev_bids.onset), 1), 'value', ev_bids.stim_type,...
            'time', num2cell(ev_bids.onset), 'duration', num2cell(ev_bids.duration));
        
        if S.replace
            D = events(D, 1, ev_spm);
        else
            D = events(D, 1, spm_cat_struct(D.events(1), ev_spm));
        end
        
        D = D.check;
        
    otherwise
        %----------------------------------------------------------------------
        fprintf('Unknown task ''%s'' to perform: Nothing done.\n',S.task);
end

% When prep is called from other functions with history, history should be
% disabled
if ~isfield(S, 'updatehistory') || S.updatehistory
    Stemp = S;
    Stemp.D = fullfile(D.path,D.fname);
    Stemp.save = 1;
    D = D.history('spm_eeg_prep', Stemp);
end

if isfield(S, 'save') && S.save
    save(D);
end


%==========================================================================
% function coreg
%==========================================================================
function M1 = coreg(fid, shape, regfid)
[junk, sel1] = spm_match_str(regfid(:, 1), fid.fid.label);
[junk, sel2] = spm_match_str(regfid(:, 2), shape.fid.label);

S = [];
S.targetfid = fid;
S.targetfid.fid.pnt = S.targetfid.fid.pnt(sel1, :);

S.sourcefid = shape;
S.sourcefid.fid.pnt = S.sourcefid.fid.pnt(sel2, :);
S.sourcefid.fid.label = S.sourcefid.fid.label(sel2);

S.targetfid.fid.label = S.sourcefid.fid.label;

S.template = 1;
S.useheadshape = 0;

M1 = spm_eeg_inv_datareg(S);
