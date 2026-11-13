function [MDP,RGB] = spm_MDP_arcade(Nr,Nc,Nd,Nb)
% Creates an MDP structure for a simple arcade game
% FORMAT [MDP,RGB] = spm_MDP_arcade(Nr,Nc,Nd,Nb)
%--------------------------------------------------------------------------
% Nr    - number of rows
% Nc    - number of columns
% Nd    - number of initial states [default: 1]
% Np    - number of bombs
%
% RGB   - display structure
%
% outcomes(1) - paddle
% outcomes(2) - ball
% outcomes(3) - target
% outcomes(4) - bomb
% outcomes(5) - background
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging

if nargin < 1, Nr  = 12; end
if nargin < 2, Nc  = 9;  end
if nargin < 3, Nd  = 2;  end
if nargin < 4, Nb  = 2;  end

% preliminaries
%--------------------------------------------------------------------------
Ns(1) = Nc;             % location of paddle
Ns(2) = Nc;             % location of ball horizontal
Ns(3) = Nr;             % location of ball verical
Ns(4) = 3;              % momentum of ball horizontal
Ns(5) = 2;              % momentum of ball verical
Ns(6) = 3;              % level of target
Ns(7) = 4;              % level of bomb
NS    = prod(Ns);       % number of states

% bomb offsets
%--------------------------------------------------------------------------
bo{1} = [ 0,0];
bo{2} = [+2,1];
bo{3} = [-4,2];
bo{4} = [-1,0];
bo{5} = [-3,1];

% likelihoods and priors
%==========================================================================
for g = 1:(Nr*Nc)
    A{g} = false(5,NS);
    A{g}(end,:) = true;
end

% proprioception
%--------------------------------------------------------------------------
P     = false(Ns(1),NS);

% transition priors
%--------------------------------------------------------------------------
B{1}  = false(NS,NS,3);

% specify likelihood and transitions for every combination of states
%--------------------------------------------------------------------------
hid   = [];              % hits
cid   = [];              % misses

x0    = fix((Nc + 1)/2); % initial location (paddle)
i0    = x0;              % initial location (ball)
j0    = 5;               % initial location
p0    = 1;               % initial momentum
q0    = 2;               % initial momentum
h0    = 3;               % initial target
u0    = 2;               % initial path
for x = 1:Ns(1)          % location of paddle
    for i = 1:Ns(2)         % location of ball horizontal
        for j = 1:Ns(3)         % location of ball verical
            for p = 1:Ns(4)         % momentum of ball horizontal
                for q = 1:Ns(5)         % momentum of ball verical
                    for h = 1:Ns(6)         % level of targets
                        for b = 1:Ns(7)         % level of bomb


                            % current state
                            %----------------------------------------------
                            reset = false;
                            s     = sub2ind(Ns,x,i,j,p,q,h,b);

                            % likelihoods
                            %==============================================

                            % for all pixels
                            %----------------------------------------------
                            for g = 1:(Nr*Nc)

                                % indices of this pixel
                                %------------------------------------------
                                ij = spm_index([Nr,Nc],g);

                                % likelihoods for target
                                %------------------------------------------
                                if ij(1) <= h
                                    A{g}(:,s) = false;
                                    A{g}(3,s) = true;
                                end

                                % likelihoods for bombs
                                %------------------------------------------
                                for n = 1:Nb
                                    ib = (i0 + bo{n}(1));
                                    jb = rem(b + bo{n}(2) - 1,Ns(7)) + 1;
                                    jb = (jb + Nr - Ns(7));
                                    if ij(1) == jb && ij(2) == ib
                                        A{g}(:,s) = false;
                                        A{g}(4,s) = true;
                                    end
                                end

                                % likelihoods for ball
                                %------------------------------------------
                                if j == ij(1) && i == ij(2)
                                    A{g}(:,s) = false;
                                    A{g}(2,s) = true;
                                end

                                % likelihoods for paddle
                                %------------------------------------------
                                if ij(1) == Nr && x == ij(2)
                                    A{g}(:,s) = false;
                                    A{g}(1,s) = true;
                                end

                            end

                            % interoception (hid and cid)
                            %==============================================

                            % hit on target row
                            %----------------------------------------------
                            if j == (h + 1) && q == 1

                                % goal
                                %------------------------------------------
                                if i > 2 && i < (Nc - 2)
                                    hid(end + 1) = s;
                                end

                            end

                            % miss
                            %----------------------------------------------
                            if j == Nr && x ~= i
                                cid(end + 1) = s;
                            end

                            % bomb
                            %----------------------------------------------
                            for n = 1:Nb
                                ib   = (i0 + bo{n}(1));
                                jb   = rem(b + bo{n}(2) - 1,Ns(7)) + 1;

                                bomb = jb == Ns(7) & x == ib;
                                bomb = bomb & ~(i == i0 & j == Nr);
                                if bomb
                                    cid(end + 1) = s;
                                end
                            end

                            % proprioception (P)
                            %==============================================
                            P(x,s) = true;

                            % transition priors
                            %==============================================
                            mp    = [-1,0,1];
                            mq    = [-1,1];
                            U     = [-1,0,1];
                            for u = 1:3

                                % next state: paddle
                                %------------------------------------------
                                xt = x + U(u);
                                if xt > Nc, xt = Nc;  end
                                if xt < 1 , xt = 1;   end

                                % next state: ball vertical
                                %------------------------------------------
                                ht = h;
                                qt = q;
                                if j <= h + 1, qt = 2; end
                                if j == Nr,    qt = 1; end

                                jt   = j + mq(qt);

                                % next state: ball horizontal
                                %------------------------------------------
                                pt = p;
                                if i == 1 && p == 1
                                    pt = 3;
                                end
                                if i == Nc && p == 3
                                    pt = 1;
                                end

                                % ball position
                                %------------------------------------------
                                it = i + mp(pt);
                                it = max(min(it,Nc),1);

                                % next state: bomb vertical
                                %------------------------------------------
                                bt = rem(b,Ns(7)) + 1;
                                
                                % miss and reset ball
                                %==========================================
                                if j == Nr                % ball at bottom

                                    if x == i

                                        % hit
                                        %----------------------------------
                                        pt = u;
                                        it = i + mp(pt);
                                        it = max(min(it,Nc),1);

                                    else

                                        % miss
                                        %----------------------------------
                                        jt = j0;
                                        pt = p0;
                                        qt = q0;

                                        % transitions
                                        %----------------------------------
                                        for d = -Nd:Nd
                                            it = max(min(i0 + d,Nc),1);
                                            st = sub2ind(Ns,xt,it,jt,pt,qt,ht,bt);
                                            B{1}(st,s,u) = true;
                                        end

                                        % reset ball
                                        %----------------------------------
                                        it = i0;

                                    end

                                end

                                % hit on target row
                                %------------------------------------------
                                if j == (h + 1) && q == 1

                                    % eliminate target row
                                    %--------------------------------------
                                    ht = max(h - 1,1);

                                end

                                % reset game
                                %==========================================

                                % hit on last row
                                %------------------------------------------
                                if j == 2 && q == 1

                                    % xt = x0;
                                    jt = j0;
                                    pt = p0;
                                    qt = q0;
                                    ht = h0;

                                    % transitions
                                    %--------------------------------------
                                    for d = -Nd:Nd
                                        it = max(min(i0 + d,Nc),1);
                                        st = sub2ind(Ns,xt,it,jt,pt,qt,ht,bt);
                                        B{1}(st,s,u) = true;
                                    end

                                    % default reset of ball
                                    %--------------------------------------
                                    it = i0;

                                end

                                % transitions
                                %------------------------------------------
                                st           = sub2ind(Ns,xt,it,jt,pt,qt,ht,bt);
                                B{1}(st,s,u) = true;
                            end

                        end
                    end
                end
            end
        end
    end
end

% control modalities (last row)
%--------------------------------------------------------------------------
for s = 1:Nc
    con(s) = sub2ind([Nr,Nc],Nr,s);
end

% ensure intened states are unqiue
%==========================================================================
hid = hid(~ismember(hid,cid));
hid = unique(hid);
cid = unique(cid);

% expose rewards, cost and action to outcomes
%==========================================================================
for g = 1:numel(A)
    id.A{g} = 1;                           % states of ball and bat
end

% hits
%--------------------------------------------------------------------------
a        = false(2,size(B{1},1));
a(1,:,:) = true;
for s = 1:size(hid,2)
    a(1,hid(1,s)) = false;
    a(2,hid(1,s)) = true;
end
A{end + 1}    = a;
id.A{end + 1} = 1;
id.reward     = numel(A);

% and misses
%--------------------------------------------------------------------------
a        = false(2,size(B{1},1));
a(1,:,:) = true;
for s = 1:size(cid,2)
    a(1,cid(1,s)) = false;
    a(2,cid(1,s)) = true;
end
A{end + 1}    = a;
id.A{end + 1} = 1;
id.contraint  = numel(A);

% and proprioception
%--------------------------------------------------------------------------
A{end + 1}    = P;
id.A{end + 1} = 1;
id.control    = numel(A);


% Enumerate the states and paths of the ensuing generative model
%==========================================================================
Ng    = numel(A);                    % number of outcome modalities
Nf    = numel(B);                    % number of hidden factors
for f = 1:Nf
    NS(f) = size(B{f},1);
    NU(f) = size(B{f},3);
end
for g = 1:Ng
    NO(g) = size(A{g},1);
end

% priors: (cost) C
%--------------------------------------------------------------------------
for g = 1:Ng
    C{g} = spm_softmax(ones(NO(g),1));
end

% This concludes the ABC of the model; namely, the likelihood mapping,
% prior transitions and preferences. Now, specify prior beliefs about
% initial states (D) and paths through those states (E)
%--------------------------------------------------------------------------
s0    = sub2ind(Ns,x0,i0,j0,p0,q0,h0);
for f = 1:Nf
    D{f} = sparse(s0,1,1,NS(f),1);    % First state
    E{f} = sparse(u0,1,1,NU(f),1);    % First path
    H{f} = [];                        % No intentions at this stage
end

% random starting positions
%--------------------------------------------------------------------------
for d = -Nd:Nd
    i = max(min(i0 + d,Nc),1);
    s = sub2ind(Ns,x0,i,j0,p0,q0,h0);
    D{1}(s) = true;
end

% specify controllable factors; here, the second factor
%--------------------------------------------------------------------------
U     = 1;                            % controllable factors

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
MDP.N = 0;                            % planning depth (1)

MDP.id = id;                          % edges


% RGB structure
%==========================================================================
% outcomes(1) - paddle
% outcomes(2) - ball
% outcomes(3) - target centre
% outcomes(4) - target periphery
% outcomes(5) - background

n     = 8;
RGB.N = [3 Nr*n Nc*n];                % image size

ball  = imread('baseball.png');
bat   = imread('paddle.png');
ufo1  = imread('ufo.png');
back  = zeros(size(ball));
ufo2  = uint8(spm_cross(sum(ufo1,3),[1 1 1]/3));

img   = {bat,ball,ufo1,ufo2,back};
for i = 1:numel(img)
    img{i} = imresize(img{i},[n,n]);
    V(:,i) = spm_vec(permute(img{i},[3,1,2]));
end

I      = cell(Nr,Nc);
[I{:}] = deal(zeros([n n]));
for i  = 1:Nr
    for j = 1:Nc

        % find indices of RGB image patches
        %-------------------------------------------------------------------
        k      = I;
        k{i,j} = k{i,j} + 1;
        k      = spm_cat(k);
        k      = spm_cross(k,ones(1,3));
        k      = permute(k,[3 1 2]);
        k      = find(k(:));

        RGB.G{i,j} = k(:);              % cell array of indices
        RGB.V{i,j} = V;                 % basis functions

    end
end


return