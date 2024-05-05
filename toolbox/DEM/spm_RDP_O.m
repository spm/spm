function MDP = spm_RDP_O(MDP,O,L)
% places outcomes in a recursive model
% FORMAT RDP = spm_RDP_O(mdp,O,[L])
% MDP - recursive MDP
% O   - outcomes
% L   - level [default: 1]
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging

% if level is not specified assume it is the first
%--------------------------------------------------------------------------
if nargin < 3, L = 1; end

% if there are multiple structures update each separately
%--------------------------------------------------------------------------
if numel(MDP) > 1
    for m = 1:size(MDP,1)                      % number of trials
        for i = 1:size(MDP,2)                  % number of agents
            mdp(m,i) = spm_RDP_O(MDP(m,i),O,L);
        end
    end
    MDP   = mdp;
    return
end


% Assume time scaling with a scale doubling
%--------------------------------------------------------------------------
str = 'MDP';
for n = 1:16
    if  eval([str '.L == L']) 
        eval([str '.S = O;'])
        return
    else
        str = [str '.MDP'];
    end
end

return