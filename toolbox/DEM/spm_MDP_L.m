function L = spm_MDP_L(P,M,U,Y)
% log-likelihood function
% FORMAT L = spm_mdp_L(P,M,U,Y)
% P    - parameter structure
% M    - generative model
% U    - inputs (observations or stimuli)
% Y    - observed responses (or choices)
%
% This auxiliary function evaluates the log likelihood of a sequence of
% choices within and between trials under and MDP model of choice behaviour
% parameterised by P.required fields of the model MR:
%
% M.G   - a function that returns a particular MDP parameterisation; i.e.,
%         MDP = M.G(P);
%__________________________________________________________________________

% place parameters in MDP
%--------------------------------------------------------------------------
if ~isstruct(P); P = spm_unvec(P,M.pE); end
mdp          = M.G(P);

% place MDP in trial structure
%--------------------------------------------------------------------------
n            = numel(U);
[MDP(1:n)]   = deal(mdp);
[MDP(1:n).o] = U{:};
[MDP(1:n).u] = Y{:};

% solve MDP and accumulate log-likelihood
%--------------------------------------------------------------------------
MDP   = spm_MDP_VB_X(MDP);
L     = 0;
for i = 1:n;
    for t = 1:size(Y{i},2)
        sub = num2cell(Y{i}(:,t));
        L   = L + log(MDP(i).P(sub{:},t));
    end
end

% cheat scaling to mimic multiple trials
%--------------------------------------------------------------------------
L = L*32;
