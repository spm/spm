function [MDP,hid,cid] = spm_MDP_pong(Nr,Nc)
% Creates an MDP structure for a simple game of Pong
% FORMAT [MDP, hid] = spm_MDP_pong(Nr,Nc)
%--------------------------------------------------------------------------
% Nr    = 6;                             % number of rows
% Nc    = 8;                             % number of columns
%
% hid   - Hidden states corresponding to hits
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging

% hid   - Hidden states corresponding to hits
%--------------------------------------------------------------------------
Ng    = Nr*Nc;                           % number of locations
Ns    = 4098;                            % maximum number of hidden states
for g = 1:Ng                             % default likelihood mapping
    A{g}        = false(5,Ns,Nc);
    A{g}(5,:,:) = true;
end
B{1}  = false(Ns,Ns);                    % transition matrices

% Likelihood and transition tensors
%--------------------------------------------------------------------------
S     = zeros(Ng,4);
i     = 2;                               % Initial location (horizontal)
j     = 2;                               % Initial location (vertical)
p     = 1;                               % Momentum (horizontal)
q     = 1;                               % Momentum (vertical)
for s = 1:Ns

    % Check whether this state has been previously visited
    %----------------------------------------------------------------------
    r      = [i,j,p,q];
    k      = ismember(S,r,'rows');
    if any(k)
        r  = find(k);
        s  = s - 1;

        B{1}(s + 1,s,1) = false;
        B{1}(r,s,1)     = true;
        B{1}(r,s,1)     = true;
        B{1}  = B{1}(1:s,1:s);
        for g = 1:Ng
            A{g} = A{g}(:,1:s,:);
        end
        break
    else
        S(s,:) = r;
    end

    % Index of latent state
    %----------------------------------------------------------------------
    n  = sub2ind([Nr,Nc],i,j);

    A{n}(1,s,:)      = true;
    A{n}(5,s,:)      = false;
    B{1}(s + 1,s,1)  = true;

    % uncomment to show orbit
    %----------------------------------------------------------------------
    % subplot(2,1,1)
    % plot(j,i,'o'), hold on
    % axis([1 Nc 1 Nr])
    % axis image, drawnow

    % Boundary conditions (switch momentum)
    %----------------------------------------------------------------------
    if ismember(i,[1,Nr]), p = -p; end
    if ismember(j,[1,Nc]), q = -q; end
    i  = i + p;
    j  = j + q;

end
Ns    = size(B{1},2);
B{1}  = B{1}(1:Ns,1:Ns);
for g = 1:Ng
    A{g} = A{g}(:,1:Ns,:);
end

% paddle (three actions)
%--------------------------------------------------------------------------
Nu    = 3;
B{2}  = false(Nc,Nc,Nu);
for u = 1:Nu
    B{2}(:,:,u) = logical(spm_speye(Nc,Nc,u - 2,2));
end

% Add paddle to likelihood mapping (observations)
%--------------------------------------------------------------------------
for s = 1:Nc
    n = sub2ind([Nr,Nc],1,s);
    A{n}(:,:,s) = false;
    A{n}(3,:,s) = true;
end

% Enumerate the states and paths of the ensuing generative model
%--------------------------------------------------------------------------
Nf    = numel(B);                    % number of hidden factors
Ng    = numel(A);                    % number of outcome modalities
for f = 1:Nf
    Ns(f) = size(B{f},1);
    Nu(f) = size(B{f},3);
end
for g = 1:Ng
    No(g) = size(A{g},1);
end


% priors: (cost) C: mild preference for not being stimulated
%--------------------------------------------------------------------------
for g = 1:Ng
    C{g} = spm_softmax([0; 1]);
end

% This concludes the ABC of the model; namely, the likelihood mapping,
% prior transitions and preferences. Now, specify prior beliefs about
% initial states (D) and paths through those states (E)
%--------------------------------------------------------------------------
for f = 1:Nf
    D{f} = sparse(1,1,1,Ns(f),1);     % First state
    E{f} = sparse(1,1,1,Nu(f),1);     % First path
    H{f} = [];                        % No intentions at this stage
end

% Ambiguity and latent states corresponding to hits
%--------------------------------------------------------------------------
hid    = [];
cid    = [];
for s1 = 1:Ns(1)
    for s2 = 1:Ns(2)

        % Render misses ambiguous
        %------------------------------------------------------------------
        if S(s1,1) == 1 && s2 ~= S(s1,2)
            cid(:,end + 1) = [s1;s2];
            for g = 1:Ng
                %%%% A{g}(:,s1,s2) = [1; 0];
            end
        end

        % And record latent states corresponding to hits
        %------------------------------------------------------------------
        if S(s1,1) == 1 && s2 == S(s1,2)
            hid(:,end + 1) = [s1;s2];
        end

    end
end

% specify controllable factors; here, the second factor
%--------------------------------------------------------------------------
U     = [0,1];                        % controllable factors

% Assemble MDP structure, with generative process
%==========================================================================
MDP.T = 8;                            % numer of moves
MDP.U = U;                            % controllable factors
MDP.A = A;                            % likelihood probabilities
MDP.B = B;                            % transition probabilities
MDP.C = C;                            % prior preferences
MDP.D = D;                            % prior over initial states
MDP.H = H;                            % prior over final states
MDP.E = E;                            % prior over initial paths
MDP.N = 0;                            % planning depth (2)


return