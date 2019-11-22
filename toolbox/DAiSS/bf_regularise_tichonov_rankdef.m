function res = bf_regularise_tichonov_rankdef(BF, S)
% Tichonov regularisation for rank deficient matrices based on the function
% contribute by Olaf Hauk
% Copyright (C) 2018 Wellcome Centre for Human Neuroimaging

% Vladimir Litvak
% $Id: bf_regularise_tichonov_rankdef.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0   
    rank = cfg_entry;
    rank.tag = 'rank';
    rank.name = 'Dimensionality';
    rank.strtype = 'n';
    rank.num = [1 1];
    rank.val = {80};
    rank.help = {'True known data rank'};
    
    lambda = cfg_entry;
    lambda.tag = 'lambda';
    lambda.name = 'Regularisation';
    lambda.strtype = 'r';
    lambda.num = [1 1];
    lambda.val = {5};
    lambda.help = {'Select the regularisation (in %)'};
    
    res      = cfg_branch;
    res.tag  = 'tichonov_rankdef';
    res.name = 'Tichonov regularisation';
    res.val  = {rank, lambda};
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

C = BF.features.(S.modality).C;
N = BF.features.(S.modality).N;

Cinv = Tikhonov_rank_def(C, S.rank, S.lambda);

features      = BF.features.(S.modality);
features.C    = C;
features.Cinv = Cinv;
features.U    = eye(size(C, 1));

res = features;