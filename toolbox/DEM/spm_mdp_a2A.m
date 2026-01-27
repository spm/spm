function MDP = spm_mdp_a2A(MDP,FIX)
% Dirichlet counts into probabilities 
% FORMAT MDP = spm_mdp_a2A(MDP,[FIX])
% MDP{n} - Cell array of MDPs
% 
% FIX.A = 0 for learning likelihoods [default: 1]
% FIX.B = 0 for learning transitions [default: 1]
%
% This auxiliary routine converts Dirichlet likelihoods and (transition)
% priors into expected probabilities. In other words, it replaces a with A
% and b with B — for a hierarchical generative model.
%
% see: spm_mdp2rdp.m
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging

% defaults
%--------------------------------------------------------------------------
if nargin < 2
    FIX.A = 1;
    FIX.B = 1;
end

% prior concentration parameters
%==========================================================================
Nm    = numel(MDP);
for n = 1:Nm

 % remove Dirichlet counts if requested
    %----------------------------------------------------------------------
    if FIX.A
        for g = 1:numel(MDP{n}.a)
            if ~isa(MDP{n}.a{g},'function_handle')
                MDP{n}.A{g} = spm_dir_norm(MDP{n}.a{g});
            end
        end
        MDP{n} = rmfield(MDP{n},'a');
    end
    if FIX.B
        for f = 1:numel(MDP{n}.b)
            if ~isa(MDP{n}.b{f},'function_handle')
                MDP{n}.B{f} = spm_dir_norm(MDP{n}.b{f});
            end
        end
        MDP{n} = rmfield(MDP{n},'b');
    end

end

return


