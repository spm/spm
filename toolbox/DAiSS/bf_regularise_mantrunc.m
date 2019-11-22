function res = bf_regularise_mantrunc(BF, S)
% User-specified dimensional reduction
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_regularise_mantrunc.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0   
    pcadim = cfg_entry;
    pcadim.tag = 'pcadim';
    pcadim.name = 'Dimensionality';
    pcadim.strtype = 'n';
    pcadim.num = [1 1];
    pcadim.val = {100};
    pcadim.help = {'User-specified number of dimensions'};
    
    res      = cfg_branch;
    res.tag  = 'mantrunc';
    res.name = 'Manual truncation';
    res.val  = {pcadim};
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

C = BF.features.(S.modality).C;
N = BF.features.(S.modality).N;

[U, alls] = svd(C);

U    = U(:,1:S.pcadim);
C    = U'*C*U; %% compact version of the covariance matrix
Cinv = pinv_plus(C);


features      = BF.features.(S.modality);
features.C    = C;
features.Cinv = Cinv;
features.U    = U;

res = features;