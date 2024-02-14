function RDP = spm_mdp2rdp(MDP,T,p)
% Converts a cell array of MDPs into a recursive MDP
% FORMAT RDP = spm_mdp2rdp(MDP,T,p)
% MDP{n} - Cell array of MDPs
%  MDP{n}.id.D - cell array of parents of D in supraordinate outcomes
%  MDP{n}.id.E - cell array of parents of E in supraordinate outcomes
% T      - path lengths [default: T = 2]
% p      - precision (rate constant) [p = 0] 
%
% RDP    - likelihood and transition tensors for generating this sequence
%  RDP.L - level or depth
%  RDP.T - time steps
%
% This auxiliary routine takes a cell array of hierarchically arranged
% MDP’s and creates a single MDP where each subordinate MDP is a field of
% the superordinate MDP. The vertical dependencies are encoded in the cell
% arrays of parents at each level of the hierarchy, where the outcomes of
% one level are the parents of the initial states and paths of the lower
% level (encoded probabilistically in the vectors D and E of the lower
% level).
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging

% Assume time scaling with a scale doubling
%--------------------------------------------------------------------------
if nargin < 2, T = 2; end
if nargin < 3, p = 0; end

% prior concentration parameters
%--------------------------------------------------------------------------
for n = 1:numel(MDP)

    % normalised likelihoods
    %----------------------------------------------------------------------
    for g = 1:numel(MDP{n}.A)
        MDP{n}.A{g} = spm_dir_norm(MDP{n}.A{g});
    end

    % normalised transitions
    %----------------------------------------------------------------------
    for f = 1:numel(MDP{n}.B)
        b     = spm_dir_norm(MDP{n}.B{f});
        B     = b;
        if p
            for t = 1:4
                for u = 1:size(B,3)
                    B(:,:,u) = B(:,:,u) + exp(-t*p)*b(:,:,u)^(t + 1);
                end
            end
        end
        MDP{n}.B{f} = B;
    end


    % remove fields
    %----------------------------------------------------------------------
    try, MDP{n} = rmfield(MDP{n},'a'); end
    try, MDP{n} = rmfield(MDP{n},'b'); end

end

% Recursively place MDP in MDP.MDP and specify path lengths
%--------------------------------------------------------------------------
SDP   = MDP{1};
SDP.T = T;
SDP.L = 1;
for n = 2:numel(MDP)
    RDP     = MDP{n};
    RDP.MDP = SDP;
    SDP     = RDP;
    SDP.T   = T;
    SDP.L   = n;
end
RDP.L = n;

return


