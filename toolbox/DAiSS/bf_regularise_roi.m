function res = bf_regularise_roi(BF, S)
% ROI regularisation
% See: Oswal et al. Optimising beamformer regions of interest analysis, Neuroimage, 2014
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Ashwini Oswal, Vladimir Litvak
% $Id: bf_regularise_roi.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0   
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
    voidef.val = {pos, radius};
    
    vois = cfg_repeat;
    vois.tag = 'vois';
    vois.name = 'ROI(s)';
    vois.num  = [1 Inf];
    vois.values = {voidef};
    vois.val  = {}; 
    vois.help = {'This makes it possible to define new VOIs when the original source space was mesh or grid.',...
        'Only the sources present in the original source space can be used at this stage'};    
    
    manual = cfg_entry;
    manual.tag = 'manual';
    manual.name = 'User-specified';
    manual.strtype = 'n';
    manual.num = [1 1];
    manual.help = {'User-specified number of dimensions'};
    
    original      = cfg_const;
    original.tag  = 'original';
    original.name = 'Keep original';
    original.help = {'Keep the original dimensionality of the ROI'};
    original.val  = {1};
    
    minka      = cfg_const;
    minka.tag  = 'minka';
    minka.name = 'Minka truncation';
    minka.help = {'Use Bayesian algorithm to determine dimensionality'};
    minka.val  = {1};
    
    dimred        = cfg_choice;
    dimred.tag    = 'dimred';
    dimred.name   = 'Dimension choice';
    dimred.values = {original, manual, minka}; 
    dimred.val    = {original};
           
    res      = cfg_branch;
    res.tag  = 'roi';
    res.name = 'ROI regularisation';
    res.val  = {vois, dimred};
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

mnipos = spm_eeg_inv_transform_points(BF.data.transforms.toMNI,  BF.sources.pos);

ind = [];
for v = 1:numel(S.voidef)
    dist = sqrt(sum((mnipos-repmat(S.voidef(v).pos, size(mnipos, 1), 1)).^2, 2));
    if S.voidef(v).radius>0
        ind = [ind find(dist<S.voidef(v).radius)];
    else
        [minval mind] = min(dist);
        if minval>20 % if there is nothing within 2cm something must be wrong
            mind = [];
            warning(['No sources were found close enough for VOI ' num2str(2)]);
        end
        ind = [ind mind];
    end        
end
            
ind = unique(ind);

if isempty(ind)
    error('No sources were found in the defined ROI');
end

BF.features.(S.modality).chanind = S.chanind;

[L, channels] = bf_fuse_lf(BF, S.modality);

LL = cat(2, L{ind});

Clf = LL*LL';

[LU, Sv, LV] = svd(Clf);

C = BF.features.(S.modality).C;
N = BF.features.(S.modality).N;

switch char(fieldnames(S.dimred))
    case 'original'
        U = LU;
    case 'manual'
        U = LU(:, 1:S.dimred.manual);
    case 'minka'          
        [M_opt,log_ev,lambda1] = spm_pca_order (LU'*C*LU, N);
        fprintf('Estimated covariance matrix order %d\n', M_opt);
 
        [U, dum] = svd(LU'*C*LU);
        U = LU*U(:, 1:M_opt);
end

C    = U'*C*U;
Cinv = pinv_plus(C);

features      = BF.features.(S.modality);
features.C    = C;
features.Cinv = Cinv;
features.U    = U;

res = features;

%%
%% plot the mean square root error here defined in the Van Veen paper as the sum of eigenvalues of components not taken over the sum of all eigenvalues
Fgraph  = spm_figure('GetWin','Graphics'); figure(Fgraph); clf
all_eig = diag(Sv).^2;
den = sum(all_eig);
all_e = [];
for n = 1:numel(all_eig)
    v = 1:n;
    num = sum(all_eig(v));
    e = (den-num)/den;
    all_e = [all_e,e];
end
semilogy(all_e,'ro-');title('error as a function of dimensionality reduction');
hold on
plot(size(U, 2)*[1 1], ylim, 'k'); 
xlabel('components');
ylabel('error defined from Van Veen paper');