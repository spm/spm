function res = bf_output_image_kurtosis(BF, S)
% Computes kurtosis image
% Copyright (C) 2019 Wellcome Centre for Human Neuroimaging

% Vladimir Litvak
% $Id: bf_output_image_kurtosis.m 7703 2019-11-22 12:06:29Z guillaume $

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
    
    method = cfg_menu;
    method.tag = 'method';
    method.name = 'Summary method';
    method.labels = {'max', 'svd'};
    method.val = {'max'};
    method.values = {'max', 'svd'};
    method.help = {'How to summarise orientations'};
    
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
    
    image_kurtosis      = cfg_branch;
    image_kurtosis.tag  = 'image_kurtosis';
    image_kurtosis.name = 'Kurtosis image';
    image_kurtosis.val  = {whatconditions, modality};
    
    res = image_kurtosis;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

D = BF.data.D;

if isfield(S.whatconditions, 'all')
    S.whatconditions.condlabel = D.condlist;
end

for i = 1:numel(S.whatconditions.condlabel)
    
    trials{i} = D.indtrial(S.whatconditions.condlabel{i}, 'GOOD');
    
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

W = BF.inverse.(S.modality).W;
nvert = numel(W);

WW = [];

spm('Pointer', 'Watch');drawnow;

spm_progress_bar('Init', nvert, ...
    sprintf('Preparing filters')); drawnow;
if nvert > 100, Ibar = floor(linspace(1, nvert,100));
else Ibar = 1:nvert; end

for i = 1:nvert
    if ~isnan(W{i})
        w    = W{i};
        
        if size(w, 1)>1
            switch S.method
                case 'max'
                    Wc          = w* BF.features.(modalities{m}).C*w';  % bf estimated source covariance matrix
                    [dum, mi]   = max(diag(Wc));
                    w = w(mi, :);
                case 'svd'
                    %% just take top pca component for now
                    Wc          = w* BF.features.(modalities{m}).C*w'; % bf estimated source covariance matrix
                    
                    [V,dum,dum]=svd(Wc);
                    w = (V(:,1)'/sqrt(size(Wc, 1)))*w;
            end
        end                
    end
    
    WW = [WW; w(:)'];
end


s2 = zeros(size(WW, 1), 1);
m4 = zeros(size(WW, 1), 1);

spm('Pointer', 'Watch');drawnow;
spm_progress_bar('Init', ntrials , 'Computing kurtosis'); drawnow;
if ntrials  > 100, Ibar = floor(linspace(1, ntrials ,100));
else Ibar = 1:ntrials; end


for i = 1:ntrials
    Y  = U'*squeeze(D(channels, ':', alltrials(i)));
    Y  = detrend(Y', 'constant')';
    
    Ys = WW*Y;
    
    s2 = s2+sum(Ys.^2, 2);
    m4 = m4+sum(Ys.^4, 2);   
     
    if ismember(i, Ibar)
        spm_progress_bar('Set', i); drawnow;
    end
end

s2 = s2/(D.nsamples*ntrials);
m4 = m4/(D.nsamples*ntrials);

k = m4 ./ s2.^2;

spm_progress_bar('Clear');

spm('Pointer', 'Arrow');drawnow;      

image(1).val   = k;
image(1).label = ['kurtosis_'  spm_file(D.fname, 'basename')];

res = image;