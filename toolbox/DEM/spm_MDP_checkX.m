function [MDP] = spm_MDP_checkX(MDP)
% MDP structure checking for XXX routines
% FORMAT [MDP] = spm_MDP_checkX(MDP,)
%
% MDP.U(1,F)            - indices of controllable paths
% MDP.T                 - number of outcomes
%
% MDP.A{g}(No(g),Ns(1),Ns(2),...)   - likelihood of outcomes given hidden states
% MDP.B{f}(Ns(f),Ns(f),Nu(f))       - transitions among hidden under MF control states
% MDP.C{g}(No(g),1)      - prior preferences over outcomes in modality g
% MDP.D{f}(Ns(f),1)      - prior probabilities over initial states
% MDP.E{f}(Nu(f),1)      - prior probabilities over initial paths
% MDP.H{f}(Ns(f),1)      - prior probabilities over final sates
%
% MDP.a{g}              - concentration parameters for A
% MDP.b{f}              - concentration parameters for B
% MDP.c{g}              - concentration parameters for C
% MDP.d{f}              - concentration parameters for D
% MDP.e{f}              - concentration parameters for E
% MDP.h{f}              - concentration parameters for H
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% deal with a sequence of trials
%==========================================================================

% if there are multiple structures check each separately
%--------------------------------------------------------------------------
if numel(MDP) > 1
    for m = 1:size(MDP,1)                      % number of trials
        for i = 1:size(MDP,2)                  % number of agents
            mdp(m,i) = spm_MDP_checkX(MDP(m,i));
        end
    end
    MDP   = mdp;
    return
end


% get (posterior or process) likelihood and priors
%--------------------------------------------------------------------------
if ~isfield(MDP,'A'), try, MDP.A = MDP.a; end, end
if ~isfield(MDP,'B'), try, MDP.B = MDP.b; end, end

% if transition priors are not specified sssume identity mappings
%--------------------------------------------------------------------------
if ~isfield(MDP,'B')
    Ns    = size(MDP.A{1},2:ndims(A));
    for f = 1:numel(Ns)
        MDP.B{f} = eye(Ns(f),Ns(f));
    end
end

% check format of likelihood and priors
%--------------------------------------------------------------------------
if ~iscell(MDP.A), MDP.A = {full(MDP.A)}; end
if ~iscell(MDP.B), MDP.B = {full(MDP.B)}; end

if isfield(MDP,'a'), if ~iscell(MDP.a), MDP.a = {full(MDP.a)}; end; end
if isfield(MDP,'b'), if ~iscell(MDP.b), MDP.b = {full(MDP.b)}; end; end


% Ensure sum to one constraints
%==========================================================================

% outcome modalities and outcomes
%--------------------------------------------------------------------------
for g = 1:numel(MDP.A)

    % ensure probabilities are normalised  : A
    %----------------------------------------------------------------------
    if isnumeric(MDP.A{g})
        MDP.A{g} = double(MDP.A{g});
        MDP.A{g} = full(spm_dir_norm(MDP.A{g}));
    end
end

% transitions, policies and states
%--------------------------------------------------------------------------
for f = 1:numel(MDP.B)

    % ensure probabilities are normalised  : B
    %----------------------------------------------------------------------
    if isnumeric(MDP.B{f})
        MDP.B{f} = double(MDP.B{f});
        MDP.B{f} = spm_dir_norm(MDP.B{f});
    end
end

% number of outcomes and latent states
%----------------------------------------------------------------------
Ng    = numel(MDP.A);                 % number of modalities
Nf    = numel(MDP.B);                 % number of factors
for g = 1:Ng
    No(g) = size(MDP.A{g},1);         % number of outcomes
end
for f = 1:Nf
    Ns(f) = size(MDP.B{f},1);         % number of hidden states
    Nu(f) = size(MDP.B{f},3);         % number of hidden paths
end

% check policy specification
%--------------------------------------------------------------------------
if ~isfield(MDP,'U')
    MDP.U = zeros(1,Nf);
end

% check preferences
%--------------------------------------------------------------------------
if ~isfield(MDP,'C')
    for g = 1:Ng
        MDP.C{g} = spm_dir_norm(ones(No(g),1));
    end
end
for g = 1:Ng
    MDP.C{g} = spm_dir_norm(MDP.C{g});
end

% check initial states
%--------------------------------------------------------------------------
if ~isfield(MDP,'D')
    for f = 1:Nf
        MDP.D{f} = spm_dir_norm(ones(Ns(f),1));
    end
else
    for f = 1:Nf
        if isempty(MDP.D{f})
            MDP.D{f} = spm_dir_norm(ones(Ns(f),1));
        end
    end
end
for f = 1:Nf
    MDP.D{f} = MDP.D{f}(:);
end

% check paths
%--------------------------------------------------------------------------
if ~isfield(MDP,'E')
    for f = 1:Nf
        MDP.E{f} = spm_dir_norm(ones(Nu(f),1));
    end
else
    for f = 1:Nf
        if isempty(MDP.E{f})
            MDP.E{f} = spm_dir_norm(ones(Nu(f),1));
        end
    end
end
for f = 1:Nf
    MDP.E{f} = MDP.E{f}(:);
end

% check domains and co-domains
%==========================================================================
% id.g      -  Indices of selected outcome modalities              
% id.A{g}   -  Indices of parents of outcomes
% id.B{f}   -  Indices of parents of states
% id.C{g}   -  Indices of parents of contraints
% id.D{f}   -  Indices of parents of initial states (next level)
% id.E{f}   -  Indices of parents of paths (next level)

% Check indices of domains and co-domains
%--------------------------------------------------------------------------
if ~isfield(MDP,'id')
    MDP.id.g  = {1:Ng};                         % Default (all modalities)
    for g = 1:Ng
        MDP.id.A{g} = 1:(ndims(MDP.A{g}) - 1);  % Default (leading factors)
        MDP.id.C{g} = [];                       % Conditional constraints
    end
end

% Check conditional constraints
%--------------------------------------------------------------------------
if ~isfield(MDP.id,'C')
    for g = 1:Ng
        MDP.id.C{g} = [];                       % Conditional constraints
    end
end

% ensure indices are row vectors id.g
%--------------------------------------------------------------------------
try
    MDP.id.g = MDP.id.g(:)';
     for g = 1:numel(MDP.id.g)
        MDP.id.g{g} = MDP.id.g{g}(:)';
    end
catch
    MDP.id.g = {1:Ng};
end

% ensure indices are row vectors id.A{g}
%--------------------------------------------------------------------------
for g = 1:Ng
    try
        MDP.id.A{g} = MDP.id.A{g}(:)';
    catch
        MDP.id.A{g} = 1:(ndims(MDP.A{g}) - 1);
    end
    if isempty(MDP.id.A{g})
        MDP.id.A{g} = 1:(ndims(MDP.A{g}) - 1);
    end
end


return