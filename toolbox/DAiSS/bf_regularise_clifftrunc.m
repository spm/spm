function res = bf_regularise_clifftrunc(BF, S)
% Regularisation based on the sudden drop-off in the covariance
% Eigenspetrum
% Copyright (C) 2021 Wellcome Centre for Human Neuroimaging

% George O'Neill
% $Id$

%--------------------------------------------------------------------------
if nargin == 0
    zthresh = cfg_entry;
    zthresh.tag = 'zthresh';
    zthresh.name = 'Z-score threshold';
    zthresh.strtype = 'r';
    zthresh.num = [1 1];
    zthresh.val = {[-1]};
    zthresh.help = {['Set the Z-score at which to detect cliff-face '...
        'in eigenspectrum, or set to below 0 to find largest jump'...
        ' (default: -1)']};
    
    res      = cfg_branch;
    res.tag  = 'clifftrunc';
    res.name = 'Eigenspectum cliff-face truncation';
    res.val  = {zthresh};
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

if ~isfield(S,'zthresh')
    S.zthresh = -1;
end

C = BF.features.(S.modality).C;
fprintf('Orignal covariance matrix condition %e\n', cond(C));

[U, alls] = svd(C);

% detect the cliff-face in the eigenspectrum, as recommended in
% Westerner et al. (NIMG 2022,
% https://doi.org/10.1016/j.neuroimage.2021.118789)
% partly using the algorithm implemented in
% https://www.fieldtriptoolbox.org/workshop/paris2019/
% handson_sourceanalysis/
% #spatial-whitening-of-the-task-data-using-the-activity-from-the-baseline

d = -diff(log10(diag(alls))); d = d./std(d);
if S.zthresh < 0
    M_opt = find(d==max(d),1,'first');
else
    M_opt = find(d==S.zthresh,1,'first');
end

fprintf('Estimated covariance matrix order %d\n', M_opt);

% Compact version of mattrix
U    = U(:,1:M_opt);
C    = U'*C*U; %% compact version of the covariance matrix
Cinv = pinv_plus(C);

fprintf('Final covariance matrix condition %e\n', cond(C));

features      = BF.features.(S.modality);
features.C    = C;
features.Cinv = Cinv;
features.U    = U;

res = features;