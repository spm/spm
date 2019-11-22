function sourcedata_robust = bf_output_sourcedata_robust(BF, S)
% Extracts source data, handling bad data segments
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_output_sourcedata_robust.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0
    method = cfg_menu;
    method.tag = 'method';
    method.name = 'Summary method';
    method.labels = {'max', 'svd', 'keep'};
    method.val = {'max'};
    method.values = {'max', 'svd', 'keep'};
    method.help = {'How to summarise sources in the ROI'};
    
    sourcedata_robust = cfg_branch;
    sourcedata_robust.tag = 'sourcedata_robust';
    sourcedata_robust.name = 'Source data (robust)';
    sourcedata_robust.val  = {method};
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;

modalities = intersect(fieldnames(BF.features), {'EEG', 'MEG', 'MEGPLANAR'});
samples    = 1:D.nsamples;


for m  = 1:numel(modalities)
    U        = BF.features.(modalities{m}).U;
    if any(any(U-eye(size(U, 1))))
        error('Projection to modes is not supported by this plug-in');
    end
    
    chanind = D.indchannel(BF.inverse.(modalities{m}).channels);
    
    bad  = badsamples(D, chanind, samples, 1);
    
    chngpnt = [1 find(any(diff(bad, [], 2)))];
    
    nsegments = length(chngpnt);
    
    id = zeros(1, length(samples));
    
    L = {};
    
    ev = [];
    for i = 1:nsegments
        goodind = ~bad(:, chngpnt(i));
        
        if i<nsegments
            id(chngpnt(i):(chngpnt(i+1)-1)) = spm_data_id(chanind(goodind));
        else
            id(chngpnt(i):end) = spm_data_id(chanind(goodind));
        end
        
        ev(i).type  = 'artefact_filterchange';
        ev(i).value = id(chngpnt(i));
        ev(i).time  = D.time(samples(chngpnt(i))) - D.time(1) + D.trialonset(1);
    end
    
    uid = unique(id);
    for i = 1:length(uid)
        id(id == uid(i)) = i;
        uid(i) = i;
    end
    
    ftdata = [];
    lbl = {};

    if isfield(BF.sources, 'voi')
        ftdata.label = BF.sources.voi.label(:);
        for v = 1:numel(ftdata.label)
            ind = find(BF.sources.voi.pos2voi == v);
            W      = cat(1, BF.inverse.(modalities{m}).W{ind});
            
            L{v}   = cat(2, BF.inverse.(modalities{m}).L{ind});
            
            switch S.method
                case 'max'
                    
                    Wc          = W* BF.features.(modalities{m}).C*W';  % bf estimated source covariance matrix
                    [dum, mi]   = max(diag(Wc));
                    
                    L{v}        = L{v}(:, mi);
                    
                case 'keep'
                    for i = 1:size(L{v}, 2)
                        lbl{end+1, 1} = [ftdata.label{v} '_' num2str(i)];
                    end
                    
            end
        end
    else
        mnipos = spm_eeg_inv_transform_points(BF.data.transforms.toMNI, BF.sources.pos);
        for i = 1:size(mnipos, 1)
            w = BF.inverse.(modalities{m}).W{i};
            if ~isnan(w)
                ftdata.label{i} = sprintf('%.2f_%.2f_%.2f', mnipos(i, :));
                L{i}            = BF.inverse.(modalities{m}).L{i};
            end
        end
        ftdata.label = ftdata.label(:);
    end
    
    if isempty(lbl)
        lbl = ftdata.label;
    end
    
    data  = nan(length(lbl), length(samples));
    C    = BF.features.(modalities{m}).C;
    
    spm('Pointer', 'Watch');drawnow;
    spm_progress_bar('Init', length(uid), ['Projecting ' modalities{m} ' data']); drawnow;
    if length(uid) > 100, Ibar = floor(linspace(1, length(uid),100));
    else Ibar = 1:length(uid); end
    
    for i = 1:length(uid)
        goodind = ~bad(:, find(id==uid(i), 1, 'first'));
        Cy    = C(goodind, goodind);
        invCy = pinv_plus(Cy);
        
        n = 1;
        for j = 1:length(ftdata.label)
            w = [];
            for k = 1:size(L{j}, 2)
                lf    = L{j}(goodind, k);
                
                w     = [w; lf'*invCy/(lf' * invCy * lf)];                               
            end
            
            data(n:(n+size(w, 1)-1), id==uid(i)) = w*D(chanind(goodind), find(id==uid(i)));
            
            n = n+size(w, 1);
        end
        if ismember(i, Ibar)
            spm_progress_bar('Set', i); drawnow;
        end
    end
    
    spm_progress_bar('Clear');
    
    ftdata.trial{1} = data;
    ftdata.time{1}  = D.time(samples);
    ftdata.label    = lbl;
    
    sourcedata_robust.(modalities{m}).ftdata = ftdata;
    sourcedata_robust.(modalities{m}).events = ev;
end
