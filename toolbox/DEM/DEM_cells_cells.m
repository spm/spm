function DEM = DEM_cells_cells

% This demo is a hierarchical extension of DEM_cells.m, where we have 16
% ensembles comprising 16 cells. Each cell has a generative model (i.e.,
% prior beliefs) about its possible local and global cell types (i.e.,
% internal, active or sensory). Given posterior beliefs about what sort of
% self it is at the local and global level, it can then predict the local
% and global intracellular signals it would expect to receive. The ensemble
% of ensembles then converges to a point attractor; where the ensemble has
% a Markov blanket and each element of the ensemble comprises a cell – that
% is itself a Markov blanket. The focus of this simulation is how the local
% level couples to the global level and vice versa. For simplicity (and
% computational expediency) we only model one ensemble at the local level
% and assume that the remaining ensembles conform to the same (local)
% dynamics. This is effectively a mean field approximation, where
% expectations of a cell in the first ensemble about its global type are
% coupled to the corresponding expectations and the ensemble level, and
% vice versa. The results of this simulation are provided in the form of a
% movie and graphs.The figure legend is included in the code below.
%
% In this example, we have used the same generative model at both levels to
% exploit the self similar hierarchical structure that emerges. However, we
% could have used different generative models at the global and local
% levels to simulate the morphogenesis of particular organelles that have a
% different form from their constituent cellular ensembles.
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_cells_cells.m 7447 2018-10-13 15:32:10Z karl $


% preliminaries
%--------------------------------------------------------------------------
clear global
rng('default')
N         = 32;                              % duration of process

% generative process and model
%==========================================================================
M(1).E.d  = 1;                               % approximation order
M(1).E.n  = 2;                               % embedding order
M(1).E.s  = 1;                               % smoothness
M(1).E.dt = 2;                               % time step

% Target morphology (for both levels)
%==========================================================================
n     = [1 6 9];                             % internal, active and sensory
r     = [0 1 3];                             % radial distance from centre
pos   = [];                                  % position
typ   = [];                                  % cell type
for j = 1:length(n)
    for i = 1:n(j)
        pos(:,end + 1) = r(j)*[sin(2*pi*i/n(j));cos(2*pi*i/n(j))];
        typ(j,end + 1) = 1;
    end
end

% Target sensation - given the canonical location and types
%--------------------------------------------------------------------------
sen   = foxhound(pos,typ);
for i = 1:3
    P.sen(:,i) = mean(sen(:,find(typ(i,:))),2);
end

% check the predicted sensations are sufficiently orthogonal for inference
%--------------------------------------------------------------------------
disp(P.sen)

% initialise expectations and action
%--------------------------------------------------------------------------
v     = typ/2 + randn(size(typ))/8;          % states (identity)                          % predicted sensations
a.pos = pos/2 + randn(size(pos))/8;          % chemotaxis
a.sec = spm_softmax(v);                      % secretion

% restriction matrix (ensuring action is strictly local)
%--------------------------------------------------------------------------
n     = sum(n);                              % total number of cells
In    = eye(n,n);
R     = spm_cat({kron(In,ones(3,2)) kron(In,ones(3,3))});
R     = kron(ones(2,1),R);


% duplicate expectations and actions
%==========================================================================
aa.l  = a;                                   % local action
aa.g  = a;                                   % global action
aa.c  = a;                                   % global action of 1st ensemble
vv.l  = v;                                   % local hidden states
vv.g  = v;                                   % global hidden states
vv.c  = v;                                   % local global hidden states
R     = kron(speye(3,3),R);                  % augment restriction

% level 1 of generative process
%--------------------------------------------------------------------------
G(1).g  = @(x,v,a,P) Gg(x,v,a,P);
G(1).v  = Gg([],[],aa,aa);
G(1).V  = exp(16);                           % precision (noise)
G(1).U  = exp(-2);                           % precision (action)
G(1).R  = R;                                 % rate matrix
G(1).pE = aa;                                % form (action)
G(1).aP = exp(-8);                           % precision (action)

% level 2; causes (action)
%--------------------------------------------------------------------------
G(2).a  = spm_vec(aa);                       % endogenous cause (action)
G(2).v  = 0;                                 % exogenous  cause
G(2).V  = exp(16);

% generative model
%==========================================================================

% level 1 of the generative model:
%--------------------------------------------------------------------------
M(1).g  = @(x,v,P) Mg([],v,P);
M(1).v  = Mg([],vv,P);
M(1).V  = exp(4);                            % precision of sensations
M(1).pE = P;

% level 2:
%--------------------------------------------------------------------------
M(2).v  = vv;
M(2).V  = exp(-4);                           % prior precision of identity

% hidden cause and prior identity expectations (and time)- none
%--------------------------------------------------------------------------
U     = zeros(spm_length(M(2).v),N);
C     = zeros(spm_length(G(2).v),N);

% assemble model structure
%--------------------------------------------------------------------------
DEM.M = M;
DEM.G = G;
DEM.C = C;
DEM.U = U;

% solve
%==========================================================================
DEM   = spm_ADEM(DEM);
spm_DEM_qU(DEM.qU,DEM.pU);

% Free energy
%--------------------------------------------------------------------------
subplot(2,2,3), plot(-DEM.J)
title('Free energy','Fontsize',16)
xlabel('time'), ylabel('Free energy')
axis square


% Graphics
%==========================================================================

% Evolution
% -------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
for t = 1:N
    
    % mmorphogenesis
    %----------------------------------------------------------------------
    subplot(2,1,1)
    v = spm_unvec(DEM.qU.a{2}(:,t),DEM.G(1).pE);
    c = spm_unvec(DEM.qU.v{2}(:,t),DEM.M(2).v);
    for i = 1:n
        for j = 1:n
            x  = v.g.pos(1,i) + v.l.pos(1,j)/8;
            y  = v.g.pos(2,i) + v.l.pos(2,j)/8;
            cl = spm_softmax(c.l(:,j),2);
            cg = spm_softmax(c.g(:,i),2);
            plot(x,y,'o','markersize',16,...
                         'LineWidth',1,...
                         'markeredgecolor',cg), hold on
            plot(x,y,'.','markersize',16,...
                         'color',cl)
        end
    end
    axis([-1 1 -1 1]*r(end))
    axis equal, box off, drawnow, hold off
    
    % movie
    %----------------------------------------------------------------------
    Mov(t) = getframe(gca);
    
    % differentiation: local, global and local expectations about global
    %----------------------------------------------------------------------
    subplot(2,3,4)
    for i = 1:n
        col = spm_softmax(c.l(:,i),2);
        plot(i,t,'.','markersize',16,'color',col), hold on, axis([0 n 0 N])
    end
    
    subplot(2,3,5)
    for i = 1:n
        col = spm_softmax(c.g(:,i),2);
        plot(i,t,'.','markersize',16,'color',col), hold on, axis([0 n 0 N])
    end
    
    subplot(2,3,6)
    for i = 1:n
        col = spm_softmax(c.c(:,i),2);
        plot(i,t,'.','markersize',16,'color',col), hold on, axis([0 n 0 N])
    end
    
end

% labels and movie
%--------------------------------------------------------------------------
subplot(2,1,1)
set(gca,'Userdata',{Mov,8})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title({'Self-organisation and differentiation';'(left click for movie)'},'FontSize',16)
xlabel('position'), ylabel('position')

subplot(2,3,4)
title('local expectations','Fontsize',16)
xlabel('cell'), ylabel('time'), box off, hold off

subplot(2,3,5)
title('global expectations','Fontsize',16)
xlabel('cell'), ylabel('time'), box off, hold off

subplot(2,3,6)
title('local about global','Fontsize',16)
xlabel('cell'), ylabel('time'), box off, hold off

% Figure legend:
% 
% This figure shows the (final) results of self organisation of an ensemble
% of cells, where each constituent of the ensemble is itself a local
% ensemble. In this example, there are 16 cells at both the global (higher)
% and local (lower) level. The upper panel shows the final disposition of
% the ensemble (of ensembles) in terms of the location of cells, and their
% differentiation (shown in colour: internal – red, active – green and
% sensory – blue). Note that there are no external states because the
% external states comprise the Markov blankets of other ensembles. Here,
% each cell is coded with two colours. The central colour corresponds to
% expectations about the type of cell in question at the local level, while
% the peripheral circle encodes expectations at the global level. The key
% thing to observe here is the emergence of a Markov blanket at both
% levels. This reflects a particular independency structure, where internal
% cells do not influence sensory (i.e. surface) cells, in virtue of their
% separation by active cells. This separation induces conditional
% independence, because of the limited range of intracellular signals (that
% fall off with a Gaussian function of Euclidean distance). The lower
% panels show the same results in a simpler format; namely, the evolution
% of subtype expectations (i.e., differentiation) at the local (left), and
% global (middle) level. The lower right panel shows the expectations of a
% single ensemble (the first) about its role at the global level. Here, the
% first ensemble is the internal state. Note the differentiation on both a
% local and global level; while local expectations about the cells’ role
% the global level converge to the same (in general) type. In these
% simulations, we used a time step of two units (of arbitrary time) and a
% second order variational filtering scheme (heuristically, this is a
% second order generalisation of extended Kalman filtering) with hidden
% states corresponding to unknown identity in terms of cell type at the
% local and global level. Please see above for more details.
%==========================================================================  


% subroutines
%==========================================================================

% Generating sensations (self-signalling and extracellular signals)
%--------------------------------------------------------------------------
function g = Gg(x,v,a,P)

a       = spm_unvec(a,P);               % action
n       = size(a.l.pos,2);              % number of cells

% local level
%--------------------------------------------------------------------------
g.l.sec = a.l.sec;                      % secretion
g.l.sen = foxhound(a.l.pos,a.l.sec);    % sensation

% global level
%--------------------------------------------------------------------------
g.g.sec = a.g.sec;                      % secretion
g.g.sen = foxhound(a.g.pos,a.g.sec);    % sensation

% global to local coupling
%--------------------------------------------------------------------------
g.c.sec = a.c.sec;                      % secretion
g.c.sen = g.g.sen(:,1)*ones(1,n);       % sensation

% local to global coupling
%--------------------------------------------------------------------------
g.g.sec(:,1) = mean(g.c.sec,2);
g.g.sen(:,1) = mean(g.c.sen,2);

return

% Generating predictionsof cell signalling and extracellular signals
%--------------------------------------------------------------------------
function g = Mg(x,v,P)

% local level
%--------------------------------------------------------------------------
p       = spm_softmax(v.l);             % expected identity
g.l.sec = p;                            % secretion
g.l.sen = P.sen*p;                      % sensation

% global level
%--------------------------------------------------------------------------
p       = spm_softmax(v.g);             % expected identity
g.g.sec = p;                            % secretion
g.g.sen = P.sen*p;                      % sensation

% local predictions of global signals
%--------------------------------------------------------------------------
p       = spm_softmax(v.c);             % expected identity
g.c.sec = p;                            % secretion
g.c.sen = P.sen*p;                      % sensation


return

% Sensed signals
%--------------------------------------------------------------------------
function sen = foxhound(x,y)
% sen = sensation
% x   = position
% y   = cell type

[m,n] = size(y);
sen   = zeros(m,n);
k     = 2;
for i = 1:n
    for j = 1:n
        
        % distance between cells (local) or ensembles (global)
        %------------------------------------------------------------------
        d = x(:,i) - x(:,j);
        d = d'*d;
        
        % signal sensed
        %------------------------------------------------------------------
        if i ~= j
           sen(:,i) = sen(:,i) + exp(-k*d)*y(:,j);
        end

    end
end

return

