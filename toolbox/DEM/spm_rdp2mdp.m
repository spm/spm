function MDP = spm_rdp2mdp(RDP,FIX)
% Converts a recursive MDP into a cell array of MDPs 
% FORMAT MDP = spm_rdp2mdp(RDP,FIX)
% MDP{n} - Cell array of MDPs
%  MDP{n}.id.D - cell array of parents of D in supraordinate outcomes
%  MDP{n}.id.E - cell array of parents of E in supraordinate outcomes
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

% Anumber of levels
%--------------------------------------------------------------------------
Nm = RDP.L;

% Check for single level models
%--------------------------------------------------------------------------
if nargin < 2
    FIX.A = 1;
    FIX.B = 1;
end

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
%--------------------------------------------------------------------------
for n = 1:Nm

    % normalise and remove fields
    %----------------------------------------------------------------------
    if FIX.A
        if isfield(MDP{n},'a')
            for g = 1:numel(MDP{n}.a)
                MDP{n}.A{g} = spm_dir_norm(MDP{n}.a{g});
            end
            MDP{n} = rmfield(MDP{n},'a');
        end
    end
    if FIX.B
        if isfield(MDP{n},'b')
            for f = 1:numel(MDP{n}.b)
                MDP{n}.B{f} = spm_dir_norm(MDP{n}.b{f});
            end
            MDP{n} = rmfield(MDP{n},'b');
        end
    end

end

return


