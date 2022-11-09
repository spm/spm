function Q = spm_gen_par(P,U)
% Generate condition specific parameters using DCM for M/EEG
% FORMAT Q = spm_gen_par(P,U)
%
% P - parameters
%   P.xc - the index of the condition of interest
% U - trial-effects
%   U.X  - between-trial effects (encodes the number of trials)
%   U.dt - time bins for within-trial effects
%
% Q - Condition specific parameters 
%__________________________________________________________________________
% Jafarian, A., Litvak, V., Cagnan, H., Friston, K.J. and Zeidman, P., 2019.
% Neurovascular coupling: insights from multi-modal dynamic causal modelling
% of fMRI and MEG. arXiv preprint arXiv:1903.07478.
%__________________________________________________________________________

% Amirhossein Jafarian
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


%--------------------------------------------------------------------------
if nargin < 2, U.X = sparse(1,0); end

% between-trial (experimental) inputs
%--------------------------------------------------------------------------
if isfield(U,'X')
    X = U.X;
else
    X = sparse(1,0);
end

if ~size(X,1)
    X = sparse(1,0);
end

% condition-specific parameters
%--------------------------------------------------------------------------
c = P.xc;
Q = spm_gen_Q(P,X(c,:));
