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
% MDP.C{G}(O,T)         - prior preferences over O outsomes in modality G
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
% MDP.alpha             - precision – action selection [16]
% MDP.beta              - precision over precision (Gamma hyperprior - [1])
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_check.m 6665 2016-01-08 21:05:13Z karl $


% deal with a sequence of trials
%==========================================================================

% if there are multiple structures check each one
%--------------------------------------------------------------------------
if length(MDP) > 1
    for i = 1:length(MDP)
        MDP(i) = spm_MDP_check(MDP(i));
    end    
    return
end

% check hierarchical MDP
%--------------------------------------------------------------------------
if isfield(MDP,'MDP')
    MDP.MDP = spm_MDP_check(MDP.MDP);
    return
end

% set up and preliminaries
%==========================================================================
try
    V = MDP.U;                      % allowable actions (1,Np)
    T = MDP.T;                      % number of transitions
catch
    try
        V = MDP.V;                  % allowable policies (T - 1,Np)
        T = size(MDP.V,1) + 1;      % number of transitions
    catch
        error('please specify sequential policies V (or U and T)')
    end
end

% numbers of transitions, policies and states
%--------------------------------------------------------------------------
Nf  = numel(MDP.B);                 % number of hidden state factors
Ng  = numel(MDP.A);                 % number of outcome factors
for f = 1:Nf
    Nu(f) = size(MDP.B{f},3);       % number of hidden controls
    Ns(f) = size(MDP.B{f},1);       % number of hidden states
end
for g = 1:Ng
    No(g) = size(MDP.A{g},1);       % number of outcomes
end
if Nf  ~= size(V,3);
    error('please ensure V(:,:,1:Nf) is conistent with MDP.B{1:Nf}')
end
if Nf  ~= numel(MDP.D);
    error('please ensure V(:,:,1:Nf) is conistent with MDP.D{1:Nf}')
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
for g = 1:Ng
    if No(g) ~= size(MDP.C{g},1);
        error(['please ensure A{' num2str(g) '} and C{' num2str(g) '} are consistent'])
    end
end

return



function A = spm_norm(A)
% normalisation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
for i = 1:size(A,2)
    for j = 1:size(A,3)
        for k = 1:size(A,4)
            for l = 1:size(A,5)
                A(:,i,j,k,l) = A(:,i,j,k,l)/sum(A(:,i,j,k,l),1);
            end
        end
    end
end

function A = spm_psi(A)
% normalisation of a probability transition rate matrix (columns)
%--------------------------------------------------------------------------
for i = 1:size(A,2)
    for j = 1:size(A,3)
        for k = 1:size(A,4)
            for l = 1:size(A,5)
                A(:,i,j,k,l) = psi(A(:,i,j,k,l)) - psi(sum(A(:,i,j,k,l)));
            end
        end
    end
end

function H = spm_ent(A)
% normalisation of a probability transition matrix (columns)
%--------------------------------------------------------------------------
for i = 1:size(A,2)
    for j = 1:size(A,3)
        for k = 1:size(A,4)
            for l = 1:size(A,5)
                H(i,j,k,l) = spm_softmax(A(:,i,j,k,l))'*A(:,i,j,k,l);
            end
        end
    end
end

function sub = spm_ind2sub(siz,ndx)
% subscripts from linear index
%--------------------------------------------------------------------------
n = numel(siz);
k = [1 cumprod(siz(1:end-1))];
for i = n:-1:1,
    vi       = rem(ndx - 1,k(i)) + 1;
    vj       = (ndx - vi)/k(i) + 1;
    sub(i,1) = vj;
    ndx      = vi;
end

