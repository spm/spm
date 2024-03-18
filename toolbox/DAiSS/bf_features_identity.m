function res = bf_features_identity(BF, S)
% Identity matrix for cases when covariance is not necessary
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2016-2023 Wellcome Centre for Human Neuroimaging


if nargin == 0
    identity      = cfg_branch;
    identity.tag  = 'identity';
    identity.name = 'Identity matrix';
    identity.val  = {};
    identity.help = {'Returns identity matrix for cases when covariance is not necessary'};
    
    res = identity;
    
    return
elseif nargin < 2
    error('Two input arguments are required');
end

features.C = eye(length(S.channels));
features.N = 1;

res = features;
