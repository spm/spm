function [Nf,Ns,Nu,Ng,No] = spm_MDP_size(mdp)
% ddimensions of MDP
% FORMAT [Nf,Ns,Nu,Ng,No] = spm_MDP_size(mdp)
% Nf  - nummber of factors
% Ns  - states per factor
% Nu  - control per factors
% Ng  - number of modalities
% No  - levels per modality
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% checks
%--------------------------------------------------------------------------
if ~isfield(mdp,'a'), mdp.a = mdp.A; end
if ~isfield(mdp,'b'), mdp.b = mdp.B; end

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
