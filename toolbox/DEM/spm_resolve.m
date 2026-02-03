function mdp = spm_resolve(mdp)
% FORMAT mdp = spm_resolve(mdp)
% mdp  - MDP structure
%
% This subroutine concerns the reduction of a likelihood tensor by moving
% posterior probability mass (Dirichlet counts) and assessing  the
% resulting possterior beliefs (about parameters) in terms of expected free
% energy which, in this instance, reduces to mutual information
%
% This subroutine compresses transition tensors of Dirichlet parameters
% summarising the joint distribution over outcomes and states. This
% compression effectively minimises expected free energy by minimising
% ambiguity; namely, the average conditional entropy over outcomes over
% states. Because this minimisation entails moving posterior probability
% mass between states, the entropy over outcomes remains the same and
% therefore minimising ambiguity corresponds to maximising the mutual
% information between states and outcomes. This can also be read as
% minimising information loss.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022 Wellcome Centre for Human Neuroimaging

% reduction operator (R): lossless compression
%==========================================================================
Ns    = size(mdp.a{1}(:,:),2);
if Ns == 1, return, end

% mutual information of likelihood
%--------------------------------------------------------------------------
M0    = spm_MDP_MI(mdp.a{g});

% initial Dirichlet counts
%--------------------------------------------------------------------------
try p = mdp.p; catch, p  = 1/32; end    

% reduce likelihood matrix
%--------------------------------------------------------------------------
for i = 1:(Ns(end) - 1)

    % evaluate reduced likelihood
    %----------------------------------------------------------------------
    M     = 0;
    for g = 1:numel(mdp.a)
        
        a{g}        = mdp.a{g};
        a{g}(:,i)   = a{g}(:,i) + a{g}(:,end) - p;
        a{g}(:,end) = p;

        % mutual information
        %------------------------------------------------------------------
        M = M + spm_MDP_MI(a{g});

    end

    % score loss of mutual information
    %----------------------------------------------------------------------
    G(i) = M - M0;

end


% ccontract likelihood matrix if there is no loss of information
%--------------------------------------------------------------------------
[G,i] = max(G); i = i(1);

if G < 0, return, end

for g = 1:numel(mdp.a)
    mdp.a{g}(:,i)     = mdp.a{g}(:,i) + mdp.a{g}(:,end) - p;
    mdp.a{g}          = mdp.a{g}(:,1:end - 1);
end
try
    mdp.b{end}(i,i,:) = mdp.b{end}(i,i,:) + mdp.b{end}(end,end,:) - p;
    mdp.b{end}        = mdp.b{end}(1:end - 1,1:end - 1,:);
end

return
