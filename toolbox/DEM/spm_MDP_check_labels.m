function [MDP] = spm_MDP_check_labels(MDP,nf,ng)
% label checking for MDP routines
% FORMAT [MDP] = spm_MDP_check_labels(MDP,nf,ng)
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
% MDP.labels             - maximum of nf and ng states and outcomes [4]
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


% deal with a sequence of trials
%==========================================================================
if nargin < 2, nf = 4; end
if nargin < 3, ng = 4; end


% if there are multiple structures check each separately
%--------------------------------------------------------------------------
if numel(MDP) > 1
    for m = 1:size(MDP,1)                      % number of trials
        for i = 1:size(MDP,2)                  % number of agents
            mdp(m,i) = spm_MDP_check_labels(MDP(m,i),nf,ng);
        end
    end
    MDP   = mdp;
    return
end


% check factors and outcome modalities have proper labels
%==========================================================================
[Nf,Ns,Nu,Ng,No] = spm_MDP_size(MDP);

for i = 1:min(Nf,nf)
    
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
for i = 1:min(Ng,ng)
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


