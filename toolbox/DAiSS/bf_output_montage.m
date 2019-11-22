function mont = bf_output_montage(BF, S)
% Generates a montage for source extraction
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_output_montage.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0
    label = cfg_entry;
    label.tag = 'label';
    label.name = 'Label';
    label.strtype = 's';
    label.help = {'Label for the VOI'};
    
    pos = cfg_entry;
    pos.tag = 'pos';
    pos.name = 'MNI coordinates';
    pos.strtype = 'r';
    pos.num = [1 3];
    pos.help = {'Locations for the VOI in MNI coordinates'};
    pos.val = {};
    
    radius = cfg_entry;
    radius.tag = 'radius';
    radius.name = 'Radius';
    radius.strtype = 'r';
    radius.num = [1 1];
    radius.val = {0};
    radius.help = {'Radius (in mm) for the VOIs (leave 0 for single point)'};
    
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
    mask.help = {'Select a mask image'};
    
    maskdef = cfg_branch;
    maskdef.tag = 'maskdef';
    maskdef.name = 'Mask VOI';
    maskdef.val  = {label, mask};
    
    vois = cfg_repeat;
    vois.tag = 'vois';
    vois.name = 'Redefine VOIs';
    vois.num  = [0 Inf];
    vois.values = {voidef, maskdef};
    vois.val  = {};
    vois.help = {'This makes it possible to define new VOIs when the original source space was mesh or grid.',...
        'Only the sources present in the original source space can be used at this stage'};
    
    method = cfg_menu;
    method.tag = 'method';
    method.name = 'Summary method';
    method.labels = {'max', 'svd', 'keep'};
    method.val = {'max'};
    method.values = {'max', 'svd', 'keep'};
    method.help = {'How to summarise sources in the ROI'};
    
    mont = cfg_branch;
    mont.tag = 'montage';
    mont.name = 'Source montage';
    mont.val  = {method, vois};
    return
elseif nargin < 2
    error('Two input arguments are required');
end
modalities = intersect(fieldnames(BF.features), {'EEG', 'MEG', 'MEGPLANAR'});

for m  = 1:numel(modalities)
    U        = BF.features.(modalities{m}).U;
    
    montage          = [];
    montage.labelorg = BF.inverse.(modalities{m}).channels;
    montage.labelorg = montage.labelorg(:);
    montage.tra      = [];
    
   chantypenew = 'LFP';
    
    if isfield(BF.inverse.(modalities{m}), 'label')
        montage.labelnew = BF.inverse.(modalities{m}).label(:);
        montage.tra = cat(1, BF.inverse.(modalities{m}).W{:})*U';
    elseif isfield(BF.sources, 'voi') || (isfield(S, 'vois') && numel(S.vois)>0)
        if isfield(BF.sources, 'voi')
            montage.labelnew = BF.sources.voi.label;
        elseif isfield(S, 'vois')
            % collect labels from voi and maks definitions
            for v = 1:numel(S.vois)
                if isfield(S.vois{v}, 'voidef')
                    montage.labelnew{v} = S.vois{v}.voidef.label;
                elseif isfield(S.vois{v}, 'maskdef')
                    montage.labelnew{v} = S.vois{v}.maskdef.label;
                end;
            end;
            mnipos = spm_eeg_inv_transform_points(BF.data.transforms.toMNI,  BF.sources.pos);
        else
            error('Don''t know what to do.');
        end
        
        lbl = {};
        
        for v = 1:numel(montage.labelnew)       % v: VOI number
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
            else
                error('Don''t know what to do.');
            end
            
            if isempty(ind) && isfield(S.vois{v}, 'voidef')
                error(['No sources were found close enough for VOI ' S.vois{v}.voidef.label]);
            elseif isempty(ind) && isfield(S.vois{v}, 'maskdef')
                error(['No sources were found close enough for VOI ' S.vois{v}.maskdef.label]);
            end
            
            W   = cat(1, BF.inverse.(modalities{m}).W{ind});
            
            switch S.method
                case 'max'
                    Wc          = W* BF.features.(modalities{m}).C*W';  % bf estimated source covariance matrix
                    [dum, mi]   = max(diag(Wc));
                    montage.tra = [montage.tra; W(mi, :)*U'];
                    
                    if isfield(BF.sources, 'voi')
                        maxpos = BF.sources.voi.pos(ind(mi), :);                        
                    else                        
                        maxpos = mnipos(ind, :);
                    end
                    disp(['Selected position for ' montage.labelnew{v} ' is [' num2str(maxpos) ']']);
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
        end;
        
        if ~isempty(lbl)
            montage.labelnew = lbl;
        end
        
    else
        mnipos = spm_eeg_inv_transform_points(BF.data.transforms.toMNI, BF.sources.pos);
        for i = 1:size(mnipos, 1)
            w = BF.inverse.(modalities{m}).W{i};
            if ~isnan(w)
                montage.labelnew{i} = sprintf('%.2f_%.2f_%.2f', mnipos(i, :));
                montage.tra = [montage.tra; w*U'];
            end
        end
        chantypenew = 'SRC';
    end
    
    montage.chantypenew = repmat({chantypenew}, length(montage.labelnew), 1);
    montage.chanunitnew = repmat({'nA*m'}, length(montage.labelnew), 1);
    
    mont.(modalities{m}) = montage;
end
