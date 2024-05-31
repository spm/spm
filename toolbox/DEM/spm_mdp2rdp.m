function RDP = spm_mdp2rdp(MDP,p,q,T,FIX)
% Converts a cell array of MDPs into a recursive MDP
% FORMAT RDP = spm_mdp2rdp(MDP,p,q,T,FIX)
% MDP{n} - Cell array of MDPs
%  MDP{n}.id.D - cell array of parents of D in supraordinate outcomes
%  MDP{n}.id.E - cell array of parents of E in supraordinate outcomes
%
% p      - likelihood concentration; e.g., p = 1/32; [default: p = 0] 
% q      - prior (transition) decay; e.g., q = 1/32;  [default: q = 0] 
% T      - path lengths [default: T = 2]
% 
% FIX.A = 0 for learning likelihoods
% FIX.B = 0 for learning transitions
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
if nargin < 5
    FIX.A = 1;
    FIX.B = 1;
end

Nm    = numel(MDP);
try p = p(1:Nm); catch, p = repmat(p(1),1,Nm); end
try q = q(1:Nm); catch, q = repmat(q(1),1,Nm); end


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
for n = 1:Nm

    % normalised likelihoods
    %----------------------------------------------------------------------
    for g = 1:numel(MDP{n}.A)
        a           = MDP{n}.A{g};
        MDP{n}.a{g} = a + p(n);
    end

    % normalised transitions
    %----------------------------------------------------------------------
    for f = 1:numel(MDP{n}.B)
        b           = MDP{n}.B{f};
        MDP{n}.b{f} = b + q(n);
    end

    % remove fields
    %----------------------------------------------------------------------
    if FIX.A
        MDP{n}.A = spm_dir_norm(MDP{n}.a);
        MDP{n}   = rmfield(MDP{n},'a');
    else
        MDP{n}   = rmfield(MDP{n},'A');
    end
    if FIX.B
        MDP{n}.B = spm_dir_norm(MDP{n}.b);
        MDP{n}   = rmfield(MDP{n},'b');
    else
        MDP{n}   = rmfield(MDP{n},'B');
    end

end

% Recursively place MDP in MDP.MDP and specify path lengths
%--------------------------------------------------------------------------
SDP   = MDP{1};
SDP.T = T;
SDP.L = 1;
for n = 2:Nm
    RDP     = MDP{n};
    RDP.MDP = SDP;
    SDP     = RDP;
    SDP.T   = T;
    SDP.L   = n;
end
RDP.L = n;

return


