function res = bf_regularise_minkatrunc(BF, S)
% Bayesian regularisation based on Minka's method
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Gareth Barnes
% $Id: bf_regularise_minkatrunc.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0   
    reduce = cfg_menu;
    reduce.tag = 'reduce';
    reduce.name = 'Reduce data dimension';
    reduce.help = {'Reduce the data to spatial modes based on Bayesian PCA'};
    reduce.labels = {'yes', 'no'};
    reduce.values = {1, 0};
    reduce.val = {1};
    
    res      = cfg_branch;
    res.tag  = 'minkatrunc';
    res.name = 'Minka truncation';
    res.val  = {reduce};
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

C = BF.features.(S.modality).C;
N = BF.features.(S.modality).N;

[U, alls] = svd(C);

[M_opt,log_ev,lambda1] = spm_pca_order (C, N);

fprintf('Estimated covariance matrix order %d\n', M_opt);


if S.reduce
    U    = U(:,1:M_opt);
    C    = U'*C*U; %% compact version of the covariance matrix
    Cinv = pinv_plus(C);
else
    C = U(:,1:M_opt)*alls(1:M_opt,1:M_opt)*U(:,1:M_opt)';
    U = eye(size(C));
    Cinv = pinv_plus(C, M_opt);
end
            
features      = BF.features.(S.modality);
features.C    = C;
features.Cinv = Cinv;
features.U    = U;

res = features;