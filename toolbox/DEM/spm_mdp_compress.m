function mdp = spm_mdp_compress(mdp)
% FORMAT mdp = spm_mdp_compress(mdp)
% mdp  - MDP structure
%   mdp.a{g} - likelihood matrices 
%   mdp.b{1} - prior transition tensor 
%
% This routine compresses the likelihood and transition tensors for a
% single factor, encoded in the a and B (Dirichlet) fields of mdp. This
% lossless compression is based upon changes in the mutual information
% (i.e., negative expected free energy) after moving Dirichlet counts among
% various rows and columns.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging

% Eliminate redundant trailing dimensions of the likelihood mapping
%--------------------------------------------------------------------------
% This subroutine compresses transition tensors of Dirichlet parameters
% summarising the joint distribution over outcomes and states. This
% compression effectively minimises expected free energy by minimising
% ambiguity; namely, the average conditional entropy over outcomes over
% states. Because this minimisation entails moving posterior probability
% mass between states, the entropy over outcomes remains the same and
% therefore minimising ambiguity corresponds to maximising the mutual
% information between states and outcomes. This can also be read as
% minimising information loss.
%--------------------------------------------------------------------------


% reduction operator (R): lossless compression
%==========================================================================
[Nf,Ns,Nu,Ng,No] = spm_MDP_size(mdp);
a     = mdp.a;
b     = mdp.b;
if Ns == 1, return, end

A     = cell(Ng,1);
for g = 1:Ng
    A{g} = spm_dir_norm(a{g});
    A{g} = A{g} > 1/32;
end

% Unique outputs
%--------------------------------------------------------------------------
[~,i,j] = unique(spm_cat(A)','rows','stable');

% Compress likelihood tensors
%--------------------------------------------------------------------------
for g = 1:Ng
    A{g}  = zeros(No(g),numel(i));
    for s = 1:Ns
        A{g}(:,j(s)) = A{g}(:,j(s)) + a{g}(:,s);
    end
end

B{1}  = zeros(numel(i),numel(i),1);
for u = 1:Nu
    for s = 1:Ns
        for r = 1:Ns

            % find empty paths
            %--------------------------------------------------------------
            if b{1}(r,s,u) > 1/32
                v  = find(~any(B{1}(:,j(s),:),1),1,'first');
                if isempty(v)
                    B{1}(j(r),j(s),end + 1) = b{1}(r,s,u);
                else
                    B{1}(j(r),j(s),v) = B{1}(j(r),j(s),v) + b{1}(r,s,u);
                end
            end
        end
    end
end


% check for lossless compression (expected free energy)
%--------------------------------------------------------------------------
if spm_MDP_MI(A) >= spm_MDP_MI(a)
    mdp.a = A;
    mdp.b = B;
end

return
