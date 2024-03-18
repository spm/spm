function RDP = spm_mdp2rdp(MDP,p,q,T)
% Converts a cell array of MDPs into a recursive MDP
% FORMAT RDP = spm_mdp2rdp(MDP,p,q,T)
% MDP{n} - Cell array of MDPs
%  MDP{n}.id.D - cell array of parents of D in supraordinate outcomes
%  MDP{n}.id.E - cell array of parents of E in supraordinate outcomes
%
% p      - likelihood concentration; e.g., p = 1/32; [default: p = 0] 
% q      - prior (transition) decay; e.g., q = 1/4;  [default: q = 0] 
% T      - path lengths [default: T = 2]
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
if nargin < 2, p = 0; end
if nargin < 3, q = 0; end
if nargin < 4, T = 2; end

% Check for single level models
%--------------------------------------------------------------------------
if numel(MDP) < 2
    RDP   = MDP{1};
    RDP.L = 1;
    RDP.T = T;
    return
end

% prior concentration parameters
%--------------------------------------------------------------------------
for n = 1:numel(MDP)

    % normalised likelihoods
    %----------------------------------------------------------------------
    for g = 1:numel(MDP{n}.A)
        a           = spm_dir_norm(MDP{n}.A{g});
        MDP{n}.A{g} = spm_dir_norm(a + p);
    end

    % normalised transitions
    %----------------------------------------------------------------------
    for f = 1:numel(MDP{n}.B)
        b     = spm_dir_norm(MDP{n}.B{f});
        B     = b;
        for t = 1:4
            for u = 1:size(B,3)
                B(:,:,u) = B(:,:,u) + (q^t)*(b(:,:,u)^(t + 1));
            end
        end
        MDP{n}.B{f} = spm_dir_norm(B);
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


