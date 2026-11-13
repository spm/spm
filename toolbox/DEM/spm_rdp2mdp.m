function MDP = spm_rdp2mdp(RDP)
% Converts a recursive MDP into a cell array of MDPs 
% FORMAT MDP = spm_rdp2mdp(RDP)
%
% RDP    - Deep (recursive) MDP (i.e., RDP.MDP.MDP...)  
%  RDP.L - level or depth
%  RDP.T - time steps
%
% MDP{n} - Cell array of MDPs
%  MDP{n}.A    - likelihood tensors 
%  MDP{n}.B    - transition tensors (logical)
%  MDP{n}.id.A - cell array of parents of A factors at the same level
%  MDP{n}.id.B - cell array of parents of B factors at the same level
%  MDP{n}.id.C - cell array of parents of C priors  at the same level
%  MDP{n}.id.D - cell array of parents of initial state at next level
%  MDP{n}.id.E - cell array of parents of initial paths at next level
%
% FIX.A = 0 for learning likelihoods
% FIX.B = 0 for learning transitions
%
% This auxiliary routine takes a single MDP where each subordinate MDP is a
% field of the superordinate MDP and creates a cell array of hierarchically
% arranged MDP’s. The vertical dependencies are encoded in the cell arrays
% of parents at each level of the hierarchy, where the outcomes of one
% level are the parents of the initial states and paths of the lower level
% (encoded probabilistically in the vectors D and E of the lower level).
%
% see: spm_mdp2rdp.m
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging

% number of levels
%--------------------------------------------------------------------------
Nm = RDP.L;

% Recursively place MDP.MDP in MDP{n}
%--------------------------------------------------------------------------
MDP   = cell(1,Nm);
for n = 1:Nm
    i      = Nm - n + 1;
    MDP{i} = RDP;
    if n < Nm
        MDP{i} = rmfield(MDP{i},'MDP');
        RDP    = RDP.MDP;
    end
end

% prior concentration parameters
%==========================================================================
for n = 1:Nm

    % normalise and remove fields
    %----------------------------------------------------------------------
    if isfield(MDP{n},'a')
        for g = 1:numel(MDP{n}.a)
            MDP{n}.A{g} = spm_dir_norm(MDP{n}.a{g});
        end
    end
    if isfield(MDP{n},'b')
        for f = 1:numel(MDP{n}.b)
            MDP{n}.B{f} = spm_dir_norm(MDP{n}.b{f});
        end
    end

end

% reinstate empty columns of transition tensors
%==========================================================================
for n = 1:Nm
    for f = 1:numel(MDP{n}.B)
        b  = MDP{n}.B{f};
        Ns = size(b,2);
        Nu = size(b,3);

        % retain precise transitions
        %------------------------------------------------------------------
        if Nu > 1 && Ns > 1 && ~isa(b,'function_handle')
            for u = 1:Nu
                for s = 1:Ns
                     b(:,s,u) = b(:,s,u) > min(b(:,s,u))*16;
                end
            end
        end

        % remove degenerate transitions
        %------------------------------------------------------------------
        if Nu > 1 && Ns > 1 && ~isa(b,'function_handle')
            for u = 1:Ns
                for s = 1:Ns
                    i = find(any(squeeze(b(u,s,:)),2),1,'first');
                    b(u,s,:) = false;
                    b(u,s,i) = true;
                end
            end
        end
        MDP{n}.B{f} = b;
    end
end


return


