function [MDP,hid,cid,con,RGB] = spm_MDP_breakout(Nr,Nc,Nd)
% Creates an MDP structure for a simple game of Breakout
% FORMAT [MDP,hid,cid,con,RGB] = spm_MDP_breakout(Nr,Nc,Nd)
%--------------------------------------------------------------------------
% Nr                     % number of rows
% Nc                     % number of columns
% Nd                     % randomness of resets
% hid   - Hidden states corresponding to hits
% cid   - Hidden states corresponding to misses
% con   - output modalities reporting control
% RGB   - display structure
%
% outcomes(1) - paddle
% outcomes(2) - ball
% outcomes(3) - target centre
% outcomes(4) - target periphery
% outcomes(5) - background
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging

if nargin < 3
    Nd  = 1;            % randomness of resets
end

% preliminaries
%--------------------------------------------------------------------------
Ns(1) = Nc;             % location of paddle
Ns(2) = Nc;             % location of ball horizontal
Ns(3) = Nr;             % location of ball verical
Ns(4) = 2;              % momentum of ball
Ns(5) = 2;              % momentum of ball
Ns(6) = 3;              % states of target
NS    = prod(Ns);       % number of states

% likelihoods and priors
%==========================================================================
for g = 1:(Nr*Nc)
    A{g} = false(5,NS);
    A{g}(end,:) = true;
end

% transition priors
%--------------------------------------------------------------------------
B{1}  = false(NS,NS,3);

% specify likelihood and transitions for every combination of states
%--------------------------------------------------------------------------
WIDTH = 1;
hid   = [];             % hits
cid   = [];             % misses

x0    = fix(Nc/2);      % initial location (paddle)
i0    = fix(Nc/2);      % initial location (ball)
j0    = 5;              % initial location
p0    = 1;              % initial momentum
q0    = 2;              % initial momentum
h0    = 3;              % initial target
u0    = 2;              % initial path
for x = 1:Ns(1)         % location of paddle
    for i = 1:Ns(2)         % location of ball horizontal
        for j = 1:Ns(3)         % location of ball verical
            for p = 1:Ns(4)         % momentum of ball
                for q = 1:Ns(5)         % momentum of ball
                    for h = 1:Ns(6)         % states of target

                        % current state
                        %--------------------------------------------------
                        s     = sub2ind(Ns,x,i,j,p,q,h);

                        % likelihoods
                        %==================================================

                        % for all pixels
                        %--------------------------------------------------
                        for g = 1:(Nr*Nc)

                            % indices of this pixel
                            %----------------------------------------------
                            ind = spm_index([Nr,Nc],g);

                            % likelihoods for paddle
                            %----------------------------------------------
                            if ind(1) == Nr
                                if x == ind(2)
                                    A{g}(1,s) = true;
                                    A{g}(5,s) = false;
                                end
                            end

                            % likelihoods for target
                            %----------------------------------------------
                            if ind(1) <= h
                                if abs(ind(2) - Nc/2) > WIDTH
                                    A{g}(3,s) = true;
                                    A{g}(5,s) = false;
                                else
                                    A{g}(4,s) = true;
                                    A{g}(5,s) = false;
                                end
                            end

                            % likelihoods for ball
                            %----------------------------------------------
                            if j == ind(1) && i == ind(2)
                                A{g}(2,s) = true;
                                A{g}(5,s) = false;
                            end

                        end

                        % transition priors
                        %==================================================
                        m     = [-1,1];
                        U     = [-1,0,1];
                        for u = 1:3


                            % next state: paddle
                            %----------------------------------------------
                            xt  = x + U(u);
                            if xt > Nc, xt = 1;  end
                            if xt < 1 , xt = Nc; end

                            % next state: ball
                            %----------------------------------------------
                            ht  = h;
                            pt  = p;
                            qt  = q;
                            if j <= h + 1, qt = 2; end
                            if j == Nr,    qt = 1; end
                            if i == 1,     pt = 2; end
                            if i == Nc,    pt = 1; end

                            it   = i + m(pt);
                            jt   = j + m(qt);

                            % miss
                            %==============================================
                            if j == Nr
                                if x ~= i
                                    
                                    jt = j0;
                                    pt = p0;
                                    qt = q0;

                                    % transitions 
                                    %--------------------------------------
                                    for r = -Nd:Nd
                                        it = i0 + r;
                                        st = sub2ind(Ns,xt,it,jt,pt,qt,ht);
                                        B{1}(st,s,u) = true;
                                    end

                                    % default reset of ball
                                    %--------------------------------------
                                    it = i0;

                                    % default reset of paddle
                                    %--------------------------------------
                                    xt = x0;

                                    % miss
                                    %--------------------------------------
                                    cid(end + 1) = s;
                                end
                            end

                            % hit
                            %----------------------------------------------
                            if j == (h + 1) && q == 1
                                ht = max(h - 1,1);
                                hid(end + 1) = s;
                            end

                            % reset game
                            %==============================================
                            if j == 2 && q == 1 && abs(i - Nc/2) > WIDTH

                                jt = j0;
                                pt = p0;
                                qt = q0;
                                ht = h0;

                                % transitions
                                %------------------------------------------
                                for r = -Nd:Nd
                                    it = i0 + r;
                                    st = sub2ind(Ns,xt,it,jt,pt,qt,ht);
                                    B{1}(st,s,u) = true;
                                end

                                % default reset of ball
                                %------------------------------------------
                                it = i0;

                                % default reset of paddle
                                %------------------------------------------
                                xt = x0;

                            end

                            % transitions
                            %----------------------------------------------
                            st           = sub2ind(Ns,xt,it,jt,pt,qt,ht);
                            B{1}(st,s,u) = true;

                        end
                    end
                end
            end
        end
    end
end

% control modalities (last row)
%----------------------------------------------------------------------
for s = 1:Nc
    con(s) = sub2ind([Nr,Nc],Nr,s);
end


% ensure intened states are unqiue
%==========================================================================
hid = unique(hid);
cid = unique(cid);

% Enumerate the states and paths of the ensuing generative model
%--------------------------------------------------------------------------
Ng    = numel(A);                    % number of outcome modalities
Nf    = numel(B);                    % number of hidden factors
for f = 1:Nf
    NS(f) = size(B{f},1);
    Nu(f) = size(B{f},3);
end
for g = 1:Ng
    No(g) = size(A{g},1);
end

% priors: (cost) C: mild preference for not being stimulated
%--------------------------------------------------------------------------
for g = 1:Ng
    C{g} = spm_softmax([0; 0; 0; 1]);
end

% This concludes the ABC of the model; namely, the likelihood mapping,
% prior transitions and preferences. Now, specify prior beliefs about
% initial states (D) and paths through those states (E)
%--------------------------------------------------------------------------
s0    = sub2ind(Ns,x0,i0,j0,p0,q0,h0);
for f = 1:Nf
    D{f} = sparse(s0,1,1,NS(f),1);    % First state
    E{f} = sparse(u0,1,1,Nu(f),1);    % First path
    H{f} = [];                        % No intentions at this stage
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


% RGB structure
%==========================================================================
% outcomes(1) - paddle
% outcomes(2) - ball
% outcomes(3) - target centre
% outcomes(4) - target periphery
% outcomes(5) - background

n     = 32;
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