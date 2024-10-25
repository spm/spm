function [O] = spm_get_sequences(MDP)
% Gets sequences generated at the deepest level of an MDP
% FORMAT [O] = spm_get_sequences(MDP)
% MDP  - Generative models (hierarchical)
%
% O    - Cell array of probabilitic outcomes for each episode
%
%--------------------------------------------------------------------------
% This auxiliary routine identifies the episodes (i.e., paths) encoded by
% states at the deepest or highest level of a hierarchical MDP.
%__________________________________________________________________________


% generate an episode encoded by the (i-th) state of the highest level
%==========================================================================
T     = 2;
O     = {};
Nm    = numel(MDP);
Ni    = size(MDP{Nm}.B{1},1);
for i = 1:Ni
    s{Nm} = i;
    for n = flip(2:Nm)
        s{n - 1} = [];
        for t = 1:size(s{n},2)

            % generate level n outcomes 
            %--------------------------------------------------------------
            No    = numel(MDP{n}.id.A);
            for g = 1:No
                j         = MDP{n}.id.A{g};
                [~,j]     = max(MDP{n}.A{g}(:,s{n}(j,t)));
                o{n}(g,t) = j;
            end

            % ensuing paths at level n - 1
            %--------------------------------------------------------------
            Ns    = numel(MDP{n - 1}.id.D);
            x     = zeros(Ns,T);
            for f = 1:Ns
                x(f,1) = o{n}(MDP{n - 1}.id.D{f},t);
                u      = o{n}(MDP{n - 1}.id.E{f},t);
                for r  = 2:T
                    [~,j]  = max(MDP{n - 1}.B{f}(:,x(f,r - 1),u));
                    x(f,r) = j;
                end

            end
            s{n - 1} = [s{n - 1} x];
        end
    end

    % final outcomes
    %----------------------------------------------------------------------
    for t = 1:size(s{1},2)
        No    = numel(MDP{1}.id.A);
        for g = 1:No
            j         = MDP{1}.id.A{g};
            O{i}{g,t} = MDP{1}.A{g}(:,s{1}(j,t));
        end
    end

end

return
