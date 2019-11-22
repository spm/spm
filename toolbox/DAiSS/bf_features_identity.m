function res = bf_features_identity(BF, S)
% Returns identity matrix for cases when covariance is not necessary
% Copyright (C) 2016 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: bf_features_identity.m 7703 2019-11-22 12:06:29Z guillaume $

%--------------------------------------------------------------------------
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