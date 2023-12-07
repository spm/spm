function [MDP] = spm_MDP_checkX(MDP)
% MDP structure checking for XXX routines
% FORMAT [MDP] = spm_MDP_checkX(MDP)
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
%
% optional:
% MDP.s(f,T)            - vector of true states - for each hidden factor
% MDP.o(g,T)            - vector of outcome     - for each outcome modality
% MDP.u(h,T - 1)        - vector of action      - for each hidden factor
%
% if c or d are not specified, they will be set to default values (of no
% preferences and uniform priors over initial steps).
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

% check sizes of Dirichlet parameterisation
%--------------------------------------------------------------------------
[Nf,Ns,Nu,Ng,No] = spm_MDP_size(MDP);

% fill in (posterior or process) likelihood and priors
%--------------------------------------------------------------------------
if ~isfield(MDP,'A'), MDP.A = MDP.a; end
if ~isfield(MDP,'B') && ~isfield(MDP,'b')
    for f = 1:numel(Ns)
       MDP.B{f} = eye(Ns(f),Ns(f));
    end
end

if ~isfield(MDP,'B'), MDP.B = MDP.b; end

% check format of likelihood and priors
%--------------------------------------------------------------------------
if ~iscell(MDP.A), MDP.A = {full(MDP.A)}; end
if ~iscell(MDP.B), MDP.B = {full(MDP.B)}; end

if isfield(MDP,'a'), if ~iscell(MDP.a), MDP.a = {full(MDP.a)}; end; end
if isfield(MDP,'b'), if ~iscell(MDP.b), MDP.b = {full(MDP.b)}; end; end


% check dimensions and orders
%==========================================================================

% number of transitions, policies and states
%--------------------------------------------------------------------------
for f = 1:Nf

    % ensure probabilities are normalised  : B
    %----------------------------------------------------------------------
    MDP.B{f} = double(MDP.B{f});
    MDP.B{f} = spm_dir_norm(MDP.B{f});

end

% numbber of outcome modalities and outcomes
%--------------------------------------------------------------------------
for g = 1:Ng

    % ensure probabilities are normalised  : A
    %----------------------------------------------------------------------
    if ~islogical(MDP.A{g})
        MDP.A{g} = double(MDP.A{g});
        MDP.A{g} = full(spm_dir_norm(MDP.A{g}));
    end

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
    MDP.C{g} = spm_dir_norm(MDP.C{g}(:));
    if No(g) ~= size(MDP.C{g},1)
        error(['please ensure A{' num2str(g) '} and C{' num2str(g) '} are consistent'])
    end
end

% check initial states
%--------------------------------------------------------------------------
if ~isfield(MDP,{'D'})
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
if Nf  ~= numel(MDP.D)
    error('please check MDP.D')
end
for f = 1:Nf
    MDP.D{f} = MDP.D{f}(:);
end

% check paths
%--------------------------------------------------------------------------
if ~isfield(MDP,{'E'})
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
if Nf  ~= numel(MDP.E)
    error('please check MDP.E')
end
for f = 1:Nf
    MDP.E{f} = MDP.E{f}(:);
end

% check initial states
%--------------------------------------------------------------------------
if ~isfield(MDP,{'H'})
    for f = 1:Nf
        MDP.H{f} = [];
    end
else
    for f = 1:Nf
        if ~isempty(MDP.H{f})
            if ~any(diff(MDP.H{f}))
                MDP.H{f} = [];
            end
        end
    end
end
if Nf  ~= numel(MDP.H)
    error('please check MDP.H')
end
for f = 1:Nf
    MDP.H{f} = MDP.H{f}(:);
end

% check initial states and internal consistency
%--------------------------------------------------------------------------
for f = 1:Nf
    if Ns(f) ~= size(MDP.D{f},1)
        error(['please ensure B{' num2str(f) '} and D{' num2str(f) '} are consistent'])
    end
end

% check paths and internal consistency
%--------------------------------------------------------------------------
for f = 1:Nf
    if Nu(f) ~= size(MDP.E{f},1)
        error(['please ensure B{' num2str(f) '} and E{' num2str(f) '} are consistent'])
    end
end

% check final states and internal consistency
%--------------------------------------------------------------------------
for f = 1:Nf
    if numel(MDP.H{f})
        if Ns(f) ~= size(MDP.H{f},1)
            error(['please ensure B{' num2str(f) '} and H{' num2str(f) '} are consistent'])
        end
    end
end


% check domains and co-domains
%==========================================================================
% id.g      -  Indices of selected output modalities              
% id.A{g}   -  Indices of likelihood factors this modality (domain)


% Check indices of domains and co-domains
%--------------------------------------------------------------------------
if ~isfield(MDP,'id')
    MDP.id.g  = {1:Ng};                         % Default (all modalities)
    for g = 1:Ng
        MDP.id.A{g} = 1:(ndims(MDP.A{g}) - 1);  % Default (leading factors)
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


% check initial states
%==========================================================================
if isfield(MDP,'s')
    if size(MDP.s,1) > numel(MDP.B)
        error('please specify an initial state MDP.s for %i factors',Nf)
    end
    f  = max(MDP.s,[],2)';
    if any(f > Ns(1:numel(f)))
        error('please ensure initial states MDP.s are consistent with MDP.B')
    end
end

% check outcomes if specified
%--------------------------------------------------------------------------
if isfield(MDP,'o')
    if numel(MDP.o)
        if size(MDP.o,1) ~= Ng
            error('please specify an outcomes MDP.o for %i modalities',Ng)
        end
        if any(max(MDP.o,[],2) > No(:))
            error('please ensure # outcomes MDP.o are consistent with MDP.A')
        end
    end
end


% check factors and outcome modalities have proper labels
%--------------------------------------------------------------------------
for i = 1:min(Nf,4)
    
    % name of factors
    %----------------------------------------------------------------------
    try
        MDP.label.factor(i);
    catch
        MDP.label.factor{i} = sprintf('factor %i',i);
    end
    
    % name of levels of each factor
    %----------------------------------------------------------------------
    for j = 1:min(Ns(i),4)
        try
            MDP.label.name{i}(j);
        catch
            MDP.label.name{i}{j} = sprintf('state %i(%i)',j,i);
        end
    end
    
    % name of actions under each factor
    %----------------------------------------------------------------------
    for j = 1:min(Nu(i),4)
        try
            MDP.label.action{i}(j);
        catch
            MDP.label.action{i}{j} = sprintf('path %i(%i)',j,i);
        end
    end
end

% name of outcomes under each modality
%--------------------------------------------------------------------------
for i = 1:min(Ng,4)
    try
        MDP.label.modality(i);
    catch
        MDP.label.modality{i} = sprintf('modality %i',i);
    end
    for j = 1:No(i)
        try
            MDP.label.outcome{i}(j);
        catch
            MDP.label.outcome{i}{j} = sprintf('outcome %i(%i)',j,i);
        end
    end
end


