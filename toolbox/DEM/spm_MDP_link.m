function [LINK,link] = spm_MDP_link(MDP)
% auxiliary function to create link (cell array)
% FORMAT [LINK,link] = spm_MDP_link(MDP)
%
% MDP.MDP  - hierarchical MDP structure
%
% LINK  - cell array of (binary) matrices linking outputs to states
% link  - (binary) matrix of non-empty links
%
% this routine assumes unique names in MDP.labels
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_link.m 7382 2018-07-25 13:58:04Z karl $

% search for matching strings identifying outcomes at the higher level with
% states at the lower level
%--------------------------------------------------------------------------
Ng    = numel(MDP.MDP.label.factor);
Nf    = numel(MDP.label.modality);
LINK  = cell(Ng,Nf);
link  = zeros(Ng,Nf);
for g = 1:Ng
    for f = 1:Nf
        state = MDP.MDP.label.name{g};
        outco = MDP.label.outcome{f};
        Ns    = numel(state);
        No    = numel(outco);
        J     = zeros(Ns,No);
        for i = 1:Ns
            for j = 1:No
                J(i,j) = strcmp(state(i),outco(j));
            end
        end
        if any(J(:))
            LINK{g,f} = J;
            link(g,f) = 1;
        end
    end
end
