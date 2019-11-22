function mont = bf_output_PLI(BF, S)
% Generates a montage for source extraction, projects the sources one by
% one and computes phase lag index on the fly, which is written out with no
% need to call bf_write. Only one VOI allowed, which can be a sphere or a
% mask image.

% Copyright (C) 2015-2017 Wellcome Trust Centre for Neuroimaging & University of
% Zurich

% Dominik R Bach
% PLI computation using code published by Gerald Cooray (2010) Karolinska 
% Institutet: EFFECT OF DIABETES MELLITUS ON HUMAN BRAIN FUNCTION 
% (https://openarchive.ki.se/xmlui/bitstream/handle/10616/40241/ram_ber%C3%A4ttelse.pdf)
% with a method based on Stam CJ, Nolte G, Daffertshofer A (Hum Brain Mapp 2007)
% 
% $Id: bf_output_PLI.m 7706 2019-11-22 16:30:29Z spm $
%--------------------------------------------------------------------------
if nargin == 0
    label = cfg_entry;
    label.tag = 'label';
    label.name = 'Label';
    label.strtype = 's';
    label.help = {'Label for the reference VOI'};
    
    pos = cfg_entry;
    pos.tag = 'pos';
    pos.name = 'MNI coordinates';
    pos.strtype = 'r';
    pos.num = [1 3];
    pos.help = {'Locations for the reference VOI in MNI coordinates'};
    pos.val = {};
    
    radius = cfg_entry;
    radius.tag = 'radius';
    radius.name = 'Radius';
    radius.strtype = 'r';
    radius.num = [1 1];
    radius.val = {0};
    radius.help = {'Radius (in mm) for the reference VOIs (leave 0 for single point)'};
    
    voidef = cfg_branch;
    voidef.tag = 'voidef';
    voidef.name = 'VOI';
    voidef.val = {label, pos, radius};
    
    mask = cfg_files;
    mask.tag = 'mask';
    mask.name = 'MNI mask';
    mask.filter = 'image';
    mask.ufilter = '.*';
    mask.num     = [1 1];
    mask.help = {'Select a reference mask image'};
    
    maskdef = cfg_branch;
    maskdef.tag = 'maskdef';
    maskdef.name = 'Mask VOI';
    maskdef.val  = {label, mask};
    
    sourcefile = cfg_files;
    sourcefile.tag = 'sourcefile';
    sourcefile.name = 'Source file';
    sourcefile.filter = 'mat';
    sourcefile.num = [1 1];
    sourcefile.help = {'Select the M/EEG data file'};
    
    sourcefiledef = cfg_branch;
    sourcefiledef.tag = 'sourcefiledef';
    sourcefiledef.name = 'Source file';
    sourcefiledef.val = {label, sourcefile};
    sourcefiledef.help = {'Select an M/EEG file with channels corresponding to reference sources. Otherwise this file must exactly correspond to the file on which the BF was computed'};
    
    vois = cfg_repeat;
    vois.tag = 'vois';
    vois.name = 'Redefine reference VOIs';
    vois.num  = [1 1];
    vois.values = {voidef, maskdef, sourcefiledef};
    vois.val  = {};
    vois.help = {'Select the source VOI from which you would like to compute the PLI'};
    
    highpass = cfg_entry;
    highpass.tag = 'highpass';
    highpass.strtype = 'r';
    highpass.name = 'Highpass cutoff';
    highpass.num = [];
    highpass.val = {0};
    highpass.help = {'Highpass filter cutoff (4th order butterworth)'};
    
    lowpass = cfg_entry;
    lowpass.tag = 'lowpass';
    lowpass.name = 'Lowpass cutoff';
    lowpass.strtype = 'r';
    lowpass.num = [];
    lowpass.val = {300};
    lowpass.help = {'Lowpass filter cutoff (4th order butterworth)'};
    
    filtermethod = cfg_branch;
    filtermethod.tag = 'filtermethod';
    filtermethod.name = 'Filter reconstructed sources';
    filtermethod.val = {highpass, lowpass};
    
    method = cfg_menu;
    method.tag = 'method';
    method.name = 'Summary method';
    method.labels = {'max', 'svd', 'keep'};
    method.val = {'svd'};
    method.values = {'max', 'svd', 'keep'};
    method.help = {'How to summarise sources in the reference ROI'};
    
    mont = cfg_branch;
    mont.tag = 'PLI';
    mont.name = 'Phase Lag Index';
    mont.val  = {method, filtermethod, vois};
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;

modalities = intersect(fieldnames(BF.features), {'EEG', 'MEG', 'MEGPLANAR'});

for m  = 1:numel(modalities)
   
    U        = BF.features.(modalities{m}).U;
    
    montage          = [];
    montage.labelorg = BF.inverse.(modalities{m}).channels;
    montage.labelorg = montage.labelorg(:);
    montage.tra      = [];
    
    chantypenew = 'LFP';
    
    v = 1; % only one VOI allowed

    if isfield(BF.inverse.(modalities{m}), 'label')
        error('Don''t know what to do. This function is supposed to work with source grids/meshes.');
    elseif isfield(BF.sources, 'voi') || (isfield(S, 'vois') && numel(S.vois)>0)
        if isfield(BF.sources, 'voi')
            montage.labelnew = BF.sources.voi.label;
        elseif isfield(S, 'vois')
            % collect labels from voi and maks definitions
            if isfield(S.vois{1}, 'voidef')
                montage.labelnew{1} = S.vois{v}.voidef.label;
            elseif isfield(S.vois{1}, 'maskdef')
                montage.labelnew{1} = S.vois{v}.maskdef.label;
            elseif isfield(S.vois{1}, 'sourcefiledef')
                montage.labelnew{1} = S.vois{1}.sourcefiledef.label;
                E = spm_eeg_load(S.vois{1}.sourcefiledef.sourcefile{1});
                if numel(E(1, :, :)) ~= numel(D(1, :, :))
                    error('Reference source and BF source files have different dimensions');
                else
                    C1 = conditions(D); C2 = conditions(E);
                    for trl = 1:ntrials(D)
                        if ~strcmpi(C1{trl}, C2{trl})
                            error('Trial labels don''t match between reference source and BF source');
                        end;
                    end;
                end;
            end;
            mnipos = spm_eeg_inv_transform_points(BF.data.transforms.toMNI,  BF.sources.pos);
        else
            error('Don''t know what to do.');
        end
        
        lbl = {};
        
        if isfield(BF.sources, 'voi')
            ind = find(BF.sources.voi.pos2voi == v);
        elseif isfield(S.vois{v}, 'maskdef') % adapted from bf_sources_voi, v115
            % find sources within mask
            V   = spm_vol(char(S.vois{v}.maskdef.mask));
            vox = spm_eeg_inv_transform_points(inv(V.mat), mnipos);
            Y   = spm_sample_vol(V, vox(:, 1),  vox(:, 2), vox(:, 3), 0);
            ind = find(~isnan(Y) & abs(Y)>0);
            
            % if no sources within mask, use sources close to mask
            if isempty(ind)
                maxdist = 5;
                fprintf('No sources within mask. Using sources up to %1.1f cm away from mask.\n', maxdist/10);
                xY.def = 'mask';
                xY.spec = V;
                [xY, XYZmm, j] = spm_ROI(xY, V);
                for k =1:size(XYZmm, 2)
                    dist = sqrt(sum((mnipos-repmat(XYZmm(:,k)', size(mnipos, 1), 1)).^2, 2));
                    [minval, vxind(k)] = min(dist);
                    if minval>maxdist,vxind(k)=NaN;end;
                end;
                vxind(isnan(vxind))=[];
                ind = unique(vxind);
            end;
            
        elseif isfield(S.vois{v}, 'voidef')
            dist = sqrt(sum((mnipos-repmat(S.vois{v}.voidef.pos, size(mnipos, 1), 1)).^2, 2));
            if S.vois{v}.voidef.radius>0
                ind = find(dist<S.vois{v}.voidef.radius);
            else
                [minval ind] = min(dist);
                if minval>20 % if there is nothing within 2cm something must be wrong
                    ind = [];
                end
            end
        elseif isfield(S.vois{v}, 'sourcefiledef')
            S.method = [];
            if nchannels(E) > 1
                for chnl = 1:nchannels(E)
                    lbl{chnl} = [montage.labelnew{1}, '_', num2str(chnl)];
                end;
            else
                lbl = montage.labelnew;
            end;
        else
            error('Don''t know what to do.');
        end
        
        if isfield(S.vois{v}, 'voidef') && isempty(ind) 
            error(['No sources were found close enough for VOI ' S.vois{v}.voidef.label]);
        elseif isfield(S.vois{v}, 'maskdef') && isempty(ind)
            error(['No sources were found close enough for VOI ' S.vois{v}.maskdef.label]);
        end
        
        if ~isfield(S.vois{v}, 'sourcefiledef')
            W   = cat(1, BF.inverse.(modalities{m}).W{ind});
            switch S.method
                case 'max'
                    Wc          = W* BF.features.(modalities{m}).C*W';  % bf estimated source covariance matrix
                    [dum, mi]   = max(diag(Wc));
                    montage.tra = [montage.tra; W(mi, :)*U'];
                case 'svd'
                    %% just take top pca component for now
                    Wc          = W* BF.features.(modalities{m}).C*W'; % bf estimated source covariance matrix
                    
                    [V,dum,dum]=svd(Wc);
                    montage.tra=[montage.tra;(V(:,1)'/sqrt(size(Wc, 1)))*W*U'];
                case 'keep'
                    montage.tra = [montage.tra; W*U'];
                    for i = 1:size(W, 1)
                        lbl{end+1, 1} = [montage.labelnew{v} '_' num2str(i)];
                    end
            end;
        else
            montage.tra = NaN(numel(lbl), numel(cat(1, BF.inverse.(modalities{m}).W{1})*U'));
        end;
    else
        error('Reference source definition is required');
    end;
    
    if ~isempty(lbl)
        montage.labelnew = lbl;
    end
    
    
    % add all other sources
    spm_progress_bar('Init', size(mnipos, 1), 'Generating source projectors');
    nsources = size(mnipos, 1);
    nRef = numel(montage.labelnew);
    
    montage.tra = [montage.tra; NaN(nsources, size(montage.tra, 2))];
    for i = 1:nsources
        w = BF.inverse.(modalities{m}).W{i};
        if ~isnan(w)
            montage.tra(nRef + i, :) = w*U';
        end
        spm_progress_bar('Set', i);
    end
    
    spm_progress_bar('Clear');
    
    montage.chantypenew = repmat({chantypenew}, length(montage.labelnew), 1);
    montage.chanunitnew = repmat({'nA*m'}, length(montage.labelnew), 1);
    
    
    % doing the actual job (this is new code, up to here was modified from
    % bf_output_montage). This works trial-by-trial, source-by-source
    % ---------------------------------------------------------------------
    
    
    spm_progress_bar('Init', ntrials(D), 'Computing PLI');
    
    PLI = NaN(nsources, nRef, ntrials(D));
    
    % set up filter (4th order bandpass Butterworth)
    [b, a] = butter(2, [S.filtermethod.highpass, S.filtermethod.lowpass]./(fsample(D)/2));
    
    for iTrl = 1:ntrials(D)
        % project and keep reference sources
        for iRef = 1:nRef
            if isfield(S.vois{1}, 'sourcefiledef')
                Dnew = filtfilt(b, a, E(iRef, :, iTrl));
            else
                Dnew = filtfilt(b, a, montage.tra(iRef, :) * D(:, :, iTrl));
            end;
            complex_ref(iRef, :) = hilbert(Dnew);
        end;
        
        % project, compute  and discard all other sources
        for iOther = 1:nsources
            Dtemp = filtfilt(b, a, montage.tra(iOther + nRef, :) *  D(:, :, iTrl));
            complex_temp = hilbert(Dtemp);
            for iRef = 1:nRef
                PLI(iOther, iRef, iTrl) = abs(mean(sign(angle(complex_temp./complex_ref(iRef, :)))));
            end;
        end;
        spm_progress_bar('Set', iTrl);
    end;
    spm_progress_bar('Clear');
    
    % average over trials
    spm_progress_bar('Init', nconditions(D), 'Averaging PLI across trials');
    
    cl = condlist(D);
    aPLI = NaN(size(PLI, 1), nRef, nconditions(D));
    
    for iCond = 1:nconditions(D)
        indx = indtrial(D, cl{iCond}, 'GOOD');
        if ~isempty(indx)
            aPLI(:, :, iCond) = mean(PLI(:, :, indx), 3);
        end;
        spm_progress_bar('Set', iCond);
    end;
    spm_progress_bar('Clear');
    
    % write into NIFTI-files
    cImage = 1;
    for iRef = 1:nRef
        for iCond = 1:nconditions(D)
            BF.output.image(cImage).val = aPLI(:, iRef, iCond);
            BF.output.image(cImage).label = ['BF_', montage.labelnew{iRef}, '_', cl{iCond}];
            cImage = cImage + 1;
        end;
    end;
    S.normalise = 'none';
    S.space = 'mni';
    bf_write_nifti(BF, S);
    
end

mont = [];
