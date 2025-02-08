function [G] = spm_rgm_group(O,dx,m)
% Compression of a (Dirichlet) probability tensor
% FORMAT [G] = spm_rgm_group(O,dx)
% O  - [No x Nt cell array] of likelhoods
% dx - upper bound on size of group
% m  - number of modalities per outcome
%
% This auxiliary routine takes a set of probability distributions assembled
% in a cell array (outcomes times the number of instances). It then
% evaluates the mutual information among outcomes. In virtue of the
% positive definiteness of this matrix, one can then use the principal
% eigenvector to identify a partition of outcomes: c.f., Spectral
% clustering under the Perron–Frobenius theorem. Outcomes that do not share
% any information with other outcomes are assigned to their own group. The
% size of each subset (i.e., group) is specified with an upper bound, which
% determines the scaling law in terms of the renormalisation group.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8454 2023-11-04 17:09:11Z karl $
%__________________________________________________________________________


% preliminaries
%--------------------------------------------------------------------------
if nargin < 2, dx = 16; end
if nargin < 3, m  = 1;  end

[No,Nt] = size(O);

% return an if no outcomes
%--------------------------------------------------------------------------
if ~No,     G = {};     return, end

% return a single group if the number of outcomes is less then dx
%--------------------------------------------------------------------------
if No < dx, G = {1:No}; return, end

% deal with multiple modalities per outcome
%--------------------------------------------------------------------------
R     = {};
for t = 1:Nt
    i = 1;
    for o = 1:m:No
        p = O{o,t};
        for r = 1:(m - 1)
            p = kron(p,O{o + r,t});
        end
        R{i,t} = p;
        i = i + 1;
    end
end

% find outcomes that change over time
%--------------------------------------------------------------------------
No    = size(R,1);
Nt    = Nt - 1;
n     = false(1,No);
for o = 1:No
    n(o) = any(diff(spm_cat(R(o,:)),[],2),'all');
end

% evaluate mutual information among outcomes
%--------------------------------------------------------------------------
MI    = zeros(No,No);
for i = 1:No
    for j = i:No
        p = 0;

        % if there is shareable information
        %------------------------------------------------------------------
        if n(i) && n(j)
            for t = 1:Nt
                p = p + R{i,t}*R{j,t}';
            end
            MI(i,j) = spm_MDP_MI(p);
            MI(j,i) = MI(i,j);
        end
    end
end

% parition into groups using principal eigenvector (Perron–Frobenius)
%--------------------------------------------------------------------------
i  = 1:No;
G  = {};
dx = fix(dx);
U  = exp(-16);
while numel(i)

    % principal eigenvector and implement bound (dx)
    %----------------------------------------------------------------------
    [e,v] = eig(MI(i,i),'nobalance');
    [~,j] = max(diag(v),[],1);
    [e,j] = sort(abs(e(:,j)),'descend');
    l     = 1:min(numel(j),dx);
    j     = j(l);

    % eliminate outcomes outwith this cluster
    %----------------------------------------------------------------------
    j(e(l) < U) = [];
    G{end + 1}  = i(j);
    i(j)        = [];

end

% deal with multiple modalities per outcome
%--------------------------------------------------------------------------
for g = 1:numel(G)
    j = (G{g} - 1)*m;
    k = [];
    for i = 1:numel(j)
       k = [k, (j(i) + (1:m))];
    end
    G{g} = k;
end



return


