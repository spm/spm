function [Nf,Ns,Nu,Ng,No] = spm_MDP_size(mdp)
% Dimensions (shape) of MDP based on generative model (a,b,...) 
% FORMAT [Nf,Ns,Nu,Ng,No] = spm_MDP_size(mdp)
% Nf  - number of factors
% Ns  - states per factor
% Nu  - control per factors
% Ng  - number of modalities
% No  - levels per modality
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging


% checks
%--------------------------------------------------------------------------
if isfield(mdp,'a')
    a = mdp.a;
else
    a = mdp.A;
end
if ~iscell(a)
    a = {a};
end
if isfield(mdp,'b')
    b = mdp.b;
else
    if isfield(mdp,'B')
        b = mdp.B;
    else
        Ns    = size(a{1},2:ndims(A));
        for f = 1:numel(Ns)
            b{f} = eye(Ns(f),Ns(f));
        end
    end
end

% sizes of factors and modilities
%--------------------------------------------------------------------------
Nf    = numel(b);                    % number of hidden factors
Ng    = numel(a);                    % number of outcome modalities
Ns    = zeros(1,Nf);                 % number of hidden states
Nu    = zeros(1,Nf);                 % number of hidden controls
No    = zeros(1,Ng);                 % number of outcomes
for f = 1:Nf
    Ns(f) = size(b{f},1);
    Nu(f) = size(b{f},3);
end
for g = 1:Ng
    No(g) = size(a{g},1);
end
