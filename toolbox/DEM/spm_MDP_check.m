function [MDP] = spm_MDP_check(MDP)
% MDP structure checking
% FORMAT [MDP] = spm_MDP_check(MDP)
%
% MDP.V(T - 1,P,F)      - P allowable policies of T moves over F factors
% or
% MDP.U(1,P,F)          - P allowable actions at each move
% MDP.T                 - number of outcomes
%
% MDP.A{G}(O,N1,...,NF) - likelihood of O outcomes given hidden states
% MDP.B{F}(NF,NF,MF)    - transitions among hidden under MF control states
% MDP.C{G}(O,T)         - prior preferences over O outcomes in modality G
% MDP.D{F}(NF,1)        - prior probabilities over initial states
%
% MDP.a{G}              - concentration parameters for A
% MDP.b{F}              - concentration parameters for B
% MDP.d{F}              - concentration parameters for D
%
% optional:
% MDP.s(F,T)            - vector of true states - for each hidden factor
% MDP.o(G,T)            - vector of outcome     - for each outcome modality
% MDP.u(F,T - 1)        - vector of action      - for each hidden factor
% MDP.w(1,T)            - vector of precisions
%
% if C or D are not specified, they will be set to default values (of no
% preferences and uniform priors over initial steps).  If there are no
% policies, it will be assumed that I = 1 and all policies (for each 
% marginal hidden state) are allowed.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_MDP_check.m 6977 2016-12-24 17:48:44Z karl $
 
 
% deal with a sequence of trials
%==========================================================================
 
% if there are multiple structures check each separately
%--------------------------------------------------------------------------
if length(MDP) > 1
    for i = 1:length(MDP)
        MDP(i) = spm_MDP_check(MDP(i));
    end    
    return
end
 
% check dimensions and orders
%==========================================================================
 
% numbers of transitions, policies and states
%--------------------------------------------------------------------------
Nf  = numel(MDP.B);                 % number of hidden state factors
Ng  = numel(MDP.A);                 % number of outcome factors
for f = 1:Nf
    Nu(f)    = size(MDP.B{f},3);    % number of hidden controls
    Ns(f)    = size(MDP.B{f},1);    % number of hidden states
    MDP.B{f} = double(MDP.B{f});
end
for g = 1:Ng
    No(g)    = size(MDP.A{g},1);    % number of outcomes
    MDP.A{g} = double(MDP.A{g});
end
 
% check policy specification
%--------------------------------------------------------------------------
try
    V = MDP.U;                      % allowable actions (1,Np)
catch
    try
        V = MDP.V;                  % allowable policies (T - 1,Np)
    catch
        
        % components allowable policies from allowable actions
        %------------------------------------------------------------------
        for f = 1:Nf
            u = 1;
            for i = 1:Nf
                if i == f
                    u = kron(1:Nu(i),u);
                else
                    u = kron(ones(1,Nu(i)),u);
                end
            end
            U(1,:,f)  = u;
        end
        MDP.U = U;
        V = MDP.U;
    end
end
 
% check policy specification
%--------------------------------------------------------------------------
if Nf  ~= size(V,3);
    error('please ensure V(:,:,1:Nf) is consistent with MDP.B{1:Nf}')
end
 
% check preferences
%--------------------------------------------------------------------------
if ~isfield(MDP,'C')
    for g = 1:Ng
        MDP.C{g} = zeros(No(g),1);
    end
end
for g = 1:Ng
    if No(g) ~= size(MDP.C{g},1);
        error(['please ensure A{' num2str(g) '} and C{' num2str(g) '} are consistent'])
    end
end
 
% check iinitial states
%--------------------------------------------------------------------------
if ~isfield(MDP,'D')
    for f = 1:Nf
        MDP.D{f} = ones(Ns(f),1);
    end
end
if Nf  ~= numel(MDP.D);
    error('please ensure V(:,:,1:Nf) is consistent with MDP.D{1:Nf}')
end
 
 
% check iinitial states and internal consistency
%--------------------------------------------------------------------------
if Nf  ~= numel(MDP.D);
    error('please ensure V(:,:,1:Nf) is consistent with MDP.D{1:Nf}')
end
for f = 1:Nf
    if Ns(f) ~= size(MDP.D{f},1);
        error(['please ensure B{' num2str(f) '} and D{' num2str(f) '} are consistent'])
    end
    if Nu(f) < max(spm_vec(V(:,:,f)));
        error(['please check V(:,:,' num2str(f) ') or U(:,:,' num2str(f) ')'])
    end
    for g = 1:Ng
        Na  = size(MDP.A{g});
        if ~all(Na(2:end) == Ns);
            error(['please ensure A{' num2str(g) '} and D{' num2str(f) '} are consistent'])
        end
    end
end
 
% check probability matrices are properly specified
%--------------------------------------------------------------------------
for f = 1:Nf
    if ~all(spm_vec(any(MDP.B{f},1)))
         error(['please check B{' num2str(f) '} for missing entries'])
    end
end
for g = 1:Ng
    if ~all(spm_vec(any(MDP.A{g},1)))
         error(['please check A{' num2str(g) '} for missing entries'])
    end
end
 
% check initial states
%--------------------------------------------------------------------------
if isfield(MDP,'s')
    if size(MDP.s,1) ~= Nf
        error('please specify an initial state MDP.s for every factor')
    end
    if any(max(MDP.s,[],2) > Ns(:))
        error('please ensure initial states MDP.s are consistent with MDP.B')
    end
end
 
% check outcomes if specified
%--------------------------------------------------------------------------
if isfield(MDP,'o')
    if size(MDP.o,1) ~= Ng
        error('please specify an outcomes MDP.o for each modality')
    end
    if any(max(MDP.o,[],2) > No(:))
        error('please ensure outcomes MDP.o are consistent with MDP.A')
    end 
end
 
% check factors and outcome modalities have proper labels
%--------------------------------------------------------------------------
for i = 1:Nf
    try
        MDP.label.factor{i};
    catch
        MDP.label.factor{i} = sprintf('factor %i',i);
    end
    for j = 1:Ns(i)
        try
            MDP.label.name{i}{j};
        catch
            MDP.label.name{i}{j} = sprintf('state %i(%i)',j,i);
        end
    end
end
for i = 1:Ng
    try
        MDP.label.modality(i);
    catch
        MDP.label.modality{i} = sprintf('modality %i',i);
    end
    for j = 1:No(i)
        try
            MDP.label.outcome{i}{j};
        catch
            MDP.label.outcome{i}{j} = sprintf('outcome %i(%i)',j,i);
        end
    end
end

% check names are specified properly
%--------------------------------------------------------------------------
if isfield(MDP,'Aname')
    if numel(MDP.Aname) ~= Ng
        error('please specify an MDP.Aname for each modality')
    end
else
    MDP.Aname = MDP.label.modality;
end
if isfield(MDP,'Bname')
    if numel(MDP.Bname) ~= Nf
        error('please specify an MDP.Bname for each factor')
    end
else
    MDP.Bname = MDP.label.factor;
end

