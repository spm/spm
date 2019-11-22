function res = bf_regularise_manual(BF, S)
% Manual specification of the regularisation parameter
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_regularise_manual.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
if nargin == 0   
    lambda = cfg_entry;
    lambda.tag = 'lambda';
    lambda.name = 'Regularisation';
    lambda.strtype = 'r';
    lambda.num = [1 1];
    lambda.val = {0};
    lambda.help = {'Select the regularisation (in %)'};
    
    res      = cfg_branch;
    res.tag  = 'manual';
    res.name = 'User-specified regularisation';
    res.val  = {lambda};
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

%%%%%%%%%
%MWW 19/11/2014
% added to be compatible with multi-class festures (see bf_features)
if isfield(S,'class'),
    features_mod      = BF.features.(S.modality).class{S.class};
else   
    features_mod      = BF.features.(S.modality);    
end;
C = features_mod.C;
%%%%%%%%%

lambda = (S.lambda/100) * trace(C)/size(C,1);
C      = C + lambda * eye(size(C));
Cinv   = pinv_plus(C);
U      = eye(size(C));

features      = BF.features.(S.modality);
features.C    = C;
features.Cinv = Cinv;
features.U    = U;

res = features;