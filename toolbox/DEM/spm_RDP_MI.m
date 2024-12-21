function [MDP] = spm_RDP_MI(MDP)
% BMR of superordinate states based upon predicted outcomes
% FORMAT [MDP] = spm_RDP_MI(MDP)
%
% This routine implements a Bayesian model reduction in which the latent
% (generalised) states that the highest level of the model are merged in a
% way that entails a minimum loss of mutual information (i.e., cross
% entropy or expected free energy) between the states of the first stream
% and predicted outcomes — now and in the future — of subsequent (trailing)
% streams. These training streams would usually include generalised actions
% (i.e., policies) and constraints reported by rewarding and costly
% outcomes at the lowest level.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8454 2023-11-04 17:09:11Z karl $
%__________________________________________________________________________


% get predictive mapping from states to outcomes in trailing streams
%==========================================================================
n     = numel(MDP);
try
    A = MDP{n}.a;
catch
    A = MDP{n}.A;
end

try
    B = MDP{n}.b{1};
catch
    B = MDP{n}.B{1};
end

% preclude ambiguous transitions
%--------------------------------------------------------------------------
Ns    = size(B,2);
Nu    = size(B,3);
for u = 1:Nu
    for s = 1:Ns
        if ~any(B(:,s,u))
            [j,i]    = max(max(squeeze(B(:,s,:)),[],2));
            B(i,s,u) = j;
        end
    end
end

A     = spm_dir_norm(A);
B     = spm_dir_norm(B);
C     = {};
for s = 2:max(MDP{1}.sB)

    % parents of factor of stream S
    %----------------------------------------------------------------------
    pD = MDP{n - 1}.id.D{MDP{n - 1}.sB == s};
    pE = MDP{n - 1}.id.E{MDP{n - 1}.sB == s};

    % outcomes predicted by first stream
    %----------------------------------------------------------------------
    ps = find(ismember([MDP{n}.id.A{:}],find(MDP{n}.sB == 1)));

    % initial states of stream S predicted by first stream
    %----------------------------------------------------------------------
    pD = intersect(ps,pD);
    pE = intersect(ps,pE);

    % map to current and subsequent (generalized) outcomes
    %----------------------------------------------------------------------
    for p = 0:1
        for u = 1:Nu
            C{end + 1,1} = A{pD}*(B(:,:,u)^p);
            C{end + 1,1} = A{pE}*(B(:,:,u)^p);
        end
    end
end


% distance matrix (i.e., normalised vectors on a hypersphere)
%--------------------------------------------------------------------------
D       = spm_information_distance(C);

% discretise and return indices of unique outcomes
%--------------------------------------------------------------------------
[~,i,j] = unique(D < sqrt(2),'rows','stable');

% restriction matrix and reduction of likelihood
%--------------------------------------------------------------------------
Ns    = numel(i);
R     = sparse(1:numel(j),j,1,numel(j),Ns);
for g = 1:numel(A)
    if MDP{n}.id.A{g} == 1
        try
            MDP{n}.a{g} = MDP{n}.a{g}*R;
        catch
            MDP{n}.A{g} = logical(MDP{n}.A{g}*R);
        end
    end
end

% transitition priors
%--------------------------------------------------------------------------
B     = zeros(Ns,Ns,Nu);
try
    for u = 1:Nu
        B(:,:,u) = R'*MDP{n}.b{1}(:,:,u)*R;
    end
    MDP{n}.b{1}  = B;
catch
    for u = 1:Nu
        B(:,:,u) = logical(R'*MDP{n}.B{1}(:,:,u)*R);
    end
    MDP{n}.B{1}  = B;
end

return


