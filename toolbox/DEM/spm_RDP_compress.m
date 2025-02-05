function [MDP] = spm_RDP_compress(MDP,R)
% BMR of superordinate states based upon predicted outcomes
% FORMAT [MDP] = spm_RDP_compress(MDP,R)
% MDP   - cell array of MDP structures
% R     - reduction operator
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


% apply reduction operators
%==========================================================================

% reduction of likelihood of children of the leading factor
%--------------------------------------------------------------------------
n     = numel(MDP);
for g = 1:numel(MDP{n}.a)
    if MDP{n}.id.A{g} == 1
        try
            MDP{n}.a{g} = MDP{n}.a{g}*R;
        catch
            MDP{n}.A{g} = spm_dir_norm(MDP{n}.A{g}*R);
        end
    end
end

% transition priors
%--------------------------------------------------------------------------
Ns  = size(R,2);
Nu  = size(MDP{n}.b{1},3);
B   = zeros(Ns,Ns,Nu);
try
    for u = 1:Nu
        B(:,:,u) = R'*MDP{n}.b{1}(:,:,u)*R;
    end
    u            = any(B,[1,2]);
    MDP{n}.b{1}  = B(:,:,u);

catch
    for u = 1:Nu
        B(:,:,u) = spm_dir_norm(R'*MDP{n}.B{1}(:,:,u)*R);
    end
    u            = any(B,[1,2]);
    MDP{n}.B{1}  = B(:,:,u);end


% propogate down hierachy
%--------------------------------------------------------------------------
for n = flip(2:numel(MDP))

    % for leading stream likehoods
    %----------------------------------------------------------------------
    for i = 1:numel(MDP{n}.a)
        g = MDP{n}.id.A{i};
        if ismember(MDP{n}.sA(g),1) && ismember(MDP{n}.sC(g),1)

            % reduce likelihood
            %--------------------------------------------------------------
            R = spm_dir_reduce(MDP{n}.a{g}');
            MDP{n}.a{g} = R'*MDP{n}.a{g};

            % and their children (D - states)
            %--------------------------------------------------------------
            for f = 1:numel(MDP{n - 1}.id.D)
                if ismember(MDP{n - 1}.id.D{f},g)

                    Ns = size(R,2);
                    Nu = size(MDP{n - 1}.b{f},3);
                    B  = zeros(Ns,Ns,Nu);
                    try
                        for u = 1:Nu
                            B(:,:,u) = R'*MDP{n - 1}.b{f}(:,:,u)*R;
                        end
                        MDP{n - 1}.b{f}  = B;
                    catch
                        for u = 1:Nu
                            B(:,:,u) = R'*MDP{n - 1}.B{f}(:,:,u)*R;
                            B(:,:,u) = spm_dir_norm(B(:,:,u));
                        end
                        MDP{n - 1}.B{f}  = B;
                    end

                    % and their likelhoods
                    %------------------------------------------------------
                    for gf = find(ismember([MDP{n - 1}.id.A{:}],f))
                        MDP{n - 1}.a{gf} = MDP{n - 1}.a{gf}*R;
                    end
                end
            end

            % or children (E - paths)
            %--------------------------------------------------------------
            for f = 1:numel(MDP{n - 1}.id.E)
                if ismember(MDP{n - 1}.id.E{f},g)

                    Nu = size(R,2);
                    Ns = size(MDP{n - 1}.b{f},2);
                    B  = zeros(Ns,Ns,Nu);
                    try
                        for s = 1:Ns
                            b        = permute(MDP{n - 1}.b{f},[1 3,2]);
                            B(:,s,:) = b(:,:,s)*R;
                        end
                        MDP{n - 1}.b{f}  = B;
                    catch
                        for s = 1:Ns
                            b        = permute(MDP{n - 1}.B{f},[1 3,2]);
                            B(:,s,:) = b(:,:,s)*R;
                        end
                        MDP{n - 1}.B{f}  = B;
                    end

                end
            end
        end
    end
end

return




