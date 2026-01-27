function MDP = spm_check_edges(RDP)
% Converts a recursive MDP into a cell array of MDPs 
% FORMAT MDP = spm_check_edges(RDP)
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
% This auxiliary routine takes a single MDP where each subordinate MDP is a
% field of the superordinate MDP and creates a cell array of hierarchically
% arranged MDP’s. The vertical dependencies are encoded in the cell arrays
% of parents at each level of the hierarchy, where the outcomes of one
% level are the parents of the initial states and paths of the lower level
% (encoded probabilistically in the vectors D and E of the lower level).
%__________________________________________________________________________
 

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

% Check edges (i.e., dependencies) are consistent with the sizes of tensors
%==========================================================================

% for each level
%--------------------------------------------------------------------------
for n = 1:Nm

    % check sizes of likelihood (A) parents
    %----------------------------------------------------------------------
    for g = 1:numel(MDP{n}.id.A)

        iA    = MDP{n}.id.A{g};
        for f = 1:numel(iA)
            try
                sA = size(MDP{n}.a{g},f + 1);
            catch
                sA = size(MDP{n}.A{g},f + 1);
            end
            try
                sB = size(MDP{n}.b{iA(f)},1);
            catch
                sB = size(MDP{n}.B{iA(f)},1);
            end
            if sA ~= sB
                error('Please check A{%i} and B{%i}',g,f)
            end
        end

    end

    % for subordinate levels
    %----------------------------------------------------------------------
    if n < Nm

        % check sizes of prior (B) parents: initial states
        %------------------------------------------------------------------
        for f = 1:numel(MDP{n}.id.D)

            iD    = MDP{n}.id.D{f};
            for g = 1:numel(iD)
                try
                    sA = size(MDP{n + 1}.a{iD(g)},1);
                catch
                    sA = size(MDP{n + 1}.A{iD(g)},1);
                end
                try
                    sB = size(MDP{n}.b{f},2);
                catch
                    sB = size(MDP{n}.B{f},2);
                end
                if sA ~= sB
                    error('Please check A{%i} and B{%i}',g,f)
                end
            end
        end

        % check sizes of prior (B) parents: initial paths
        %------------------------------------------------------------------
        for f = 1:numel(MDP{n}.id.E)

            iE    = MDP{n}.id.E{f};
            for g = 1:numel(iE)
                try
                    sA = size(MDP{n + 1}.a{iE(g)},1);
                catch
                    sA = size(MDP{n + 1}.A{iE(g)},1);
                end
                try
                    sB = size(MDP{n}.b{f},3);
                catch
                    sB = size(MDP{n}.B{f},3);
                end
                if sA ~= sB
                    error('Please check A{%i} and B{%i}',g,f)
                end
            end
        end
    end
end

return


 
