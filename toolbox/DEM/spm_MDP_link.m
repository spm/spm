function [link,LINK] = spm_MDP_link(MDP)
% auxiliary function to create link (cell array)
% FORMAT [link,LINK] = spm_MDP_link(MDP)
%
% MDP.MDP  - hierarchical MDP structure
%
% link  - matrix of links
% LINK  - cell array of matrices linking outputs to states
%
% this routine assumes unique names in MDP.labels
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_link.m 7329 2018-06-10 21:12:02Z karl $

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
