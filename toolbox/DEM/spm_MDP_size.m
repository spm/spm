function [Nf,Ns,Nu,Ng,No] = spm_MDP_size(mdp)
% Dimensions of MDP
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
if ~isfield(mdp,'a'), mdp.a = mdp.A; end
if ~isfield(mdp,'b')
    Ns = size(mdp.a{1},2:16);
    Ns = Ns(Ns > 1);
    for f = 1:numel(Ns)
       mdp.b{f} = eye(Ns(f),Ns(f));
    end
end

% sizes of factors and modilities
%--------------------------------------------------------------------------
Nf    = numel(mdp.b);                    % number of hidden factors
Ng    = numel(mdp.a);                    % number of outcome modalities
Ns    = zeros(1,Nf);
Nu    = zeros(1,Nf);
No    = zeros(1,Ng);
for f = 1:Nf
    Ns(f) = size(mdp.b{f},1);            % number of hidden states
    Nu(f) = size(mdp.b{f},3);            % number of hidden controls
end
for g = 1:Ng
    No(g) = size(mdp.a{g},1);            % number of outcomes
end
