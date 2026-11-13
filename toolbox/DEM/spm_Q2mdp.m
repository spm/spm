function [mdp] = spm_Q2mdp(Q,n)
% Dimensions (shape) of MDP based on generative model (a,b,...) 
% FORMAT [mdp] = spm_Q2mdp(Q,n)
% s: states
% u: paths
% X: states (posterior)
% Y: outcomes (predictive posterior)
% O: outcomes (likelihood)
% o: outcomes
% j: outcomes (parents)
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging


% checks
%--------------------------------------------------------------------------
if nargin < 2, n = 1; end

mdp.s = Q.s{n};
mdp.u = Q.u{n};
for f = 1:size(Q.X{n},1)
    mdp.X{f,1} = spm_cat(Q.X{n}(f,:));
end
mdp.Y = Q.Y{n};
mdp.O = Q.O{n};
mdp.o = Q.o{n};
mdp.j = Q.j{n};

% sizes of factors and modilities
%--------------------------------------------------------------------------
return