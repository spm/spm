function FEP_MB_demo
% This demonstration routine simulates the emergence of life - as defined
% in terms of active inference - using a synthetic primordial soup. The key
% aspect of this dynamics is that there is a separation between dynamical
% states and structural states; where the dynamical states of the
% microsystem are equipped with a Lorentz attractor and the structural
% states correspond to position and velocity. The flow of structural
% states conforms to classical Newtonian mechanics. Crucially, the physical
% motion of each microsystem is coupled to its internal dynamics and vice
% versa; where the coupling among dynamics rests upon short range
% (electrochemical) forces. This means that the dependencies among the
% dynamics of each microsystem dependent on their positions. This induces a
% dependency of the systems structural integrity on its internal dynamics -
% which leads to biological self-organisation. This biological self-
% organisation is illustrated in terms of the following:
%
% i) the existence of a Markov blanket that separates internal and external
% states, where the internal states are associated with a system that
% engages in active or embodied inference.
%
% ii) emergent inference is demonstrated by showing that the internal
% states can predict the extent states, despite their separation by the
% Markov blanket.
%
% iii) this inference (encoded by the internal dynamics) is necessary to
% maintain structural integrity, as illustrated by simulated lesion
% experiments, in which the influence of various states are quenched.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: FEP_MB_demo.m 7032 2017-02-27 19:47:34Z karl $


% default settings (GRAPHICS sets movies)
%--------------------------------------------------------------------------
% rng('default')

% create an adjacency matrix or Jacobian based upon a lattice
%==========================================================================

% within blanket coupling (intrinsic)
%--------------------------------------------------------------------------
n     = 2;
m     = 3;
Jii   = spm_cat({[], eye(n,n),  spm_speye(n,m);
                 randn(n,n)/8, [], [];
                 [], spm_speye(m,n), (randn(m,m)/8 - eye(m,m))});


% between blanket coupling (extrinsic)
%--------------------------------------------------------------------------
Jij   = spm_cat({[], zeros(n,n),  zeros(n,m);
                 randn(n,n),  [], [];
                 [], zeros(m,n), zeros(m,m)});


% an ensemble of blankets
%--------------------------------------------------------------------------
D     = 2;
N     = 13;                      % size of lattice
[I,J] = ndgrid(1:N,1:N);        % locations
for i = 1:numel(I)
    for j = 1:numel(J)
        d = sqrt( (I(i) - I(j))^2 + (J(i) - J(j))^2 );
        if i == j
            A{i,j} = Jii;
        elseif d < D
            A{i,j} = Jij;
        else
            A{i,j} = zeros(size(Jij));
        end
    end
end
clear J
J{1}    = spm_cat(A);
z{1}    = num2cell(1:length(J{1}));


% hierarchal decomposition
%==========================================================================
N     = 3;                       % number of hierarchies
n     = [4 3 2 1];               % number of eigensolutions
m     = [3 2 1 1];               % number of internal states
for i = 1:N
    
    % Markov blanket partition
    %----------------------------------------------------------------------
    spm_figure('getwin','Markov partition');
    subplot(N,1,i)
    
    x{i} = spm_Markov_blanket(J{i},z{i},m(i));
    j    = spm_vec(x{i}');
    k    = spm_vec(x{i});
    
    % dimension reduction (eliminating internal states)
    %----------------------------------------------------------------------
    if i < N
        [J{i + 1},z{i + 1}] = spm_A_reduce(J{i},x{i},n(i));
    end

    % parameters
    %----------------------------------------------------------------------
    spm_figure('getwin','Jacobians');
    
    subplot(N,2,(i - 1)*2 + 1),imagesc(abs(J{i}(k,k))),axis square
    subplot(N,2,(i - 1)*2 + 2),imagesc(abs(J{i}(j,j))),axis square
end



% Markov blanket - parents, children, and parents of children
%==========================================================================
function [J,z] = spm_A_reduce(J,x,n)
% reduction of Markovian partition
% J  - Jacobian
% x  - {3 x N} indices of Markovian partition
% n  - number of eigenvectors to retain [default: 2]
%
% z  - {1 x N} indices of partition for the next level
%__________________________________________________________________________

% preliminaries
%--------------------------------------------------------------------------
nx    = size(x,2);                % number of partitions
if nargin < 3
    n = 2;                        % number of generalised eigenvectors
end

% reduction
%----------------------------------------------------------------------
for i = 1:nx
    Jii   = full(J(spm_vec(x(1:2,i)),spm_vec(x(1:2,i))));
    [e,s] = eig(Jii);
    [s,j] = sort(real(diag(s)),'descend');
    v{i}  = e(:,j(1:n));
    u{i}  = pinv(v{i});
end
for i = 1:nx
    for j = 1:nx
        Jij    = full(J(spm_vec(x(1:2,i)),spm_vec(x(1:2,j))));
        A{i,j} = u{i}*Jij*v{j};
    end
    z{i} = (1:n) + (i - 1)*n;
end
J     = spm_cat(A);



function [x] = spm_Markov_blanket(J,z,m)
% Markovian partition
% J  - Jacobian
% z  - {1 x N} indices of partition
% m  - number of internal states [default: 3]
%
% x  - {3 x N} indices of Markovian partition
%__________________________________________________________________________

% preliminaries
%--------------------------------------------------------------------------
nz    = length(z);                % number of partitions
if nargin < 3
    m = 3;                        % maximum size of internal states
end

% Adjacency matrix (over z)
%--------------------------------------------------------------------------
for i = 1:nz
    for j = 1:nz
        Lij    = J(z{i},z{j});
        L(i,j) = any(abs(Lij(:)) > 1e-18);
    end
end
L     = double(L);


% internal states (defined by graph Laplacian)
%--------------------------------------------------------------------------
G     = L - diag(diag(L));
G     = G - diag(sum(G));
G     = expm(G);

% get principal dimensions of scaling space (X)
%--------------------------------------------------------------------------
GRAPHICS = 1;
if GRAPHICS
    
    [e v] = eig(G,'nobalance');
    [v,j] = sort(real(diag(v)),'descend');
    X     = real(e(:,j(2:4)));
    
    % get first two dimensions of scaling space (X)
    %----------------------------------------------------------------------
    plot3(X(:,1),X(:,2),X(:,3),'.c','MarkerSize',24), hold on
    
end

% get Markov blanket and divide into sensory and active states
%--------------------------------------------------------------------------
B     = L + L' + L*L';
B     = B - diag(diag(B));
nn    = zeros(nz,1);

% recursive partition
%--------------------------------------------------------------------------
for i = 1:128
    
    % internal states (defined by graph Laplacian)
    %----------------------------------------------------------------------
    [g,j] = max(max(G)'.*~(B*nn));
    if g == 0
        break
    end
    g     = G(:,j);
    [g,j] = sort(g,'descend');   
    j     = j(1:m);                                   % internal states

    jj    = sparse(j,1,1,size(L,1),1);                % internal states
    bb    = B*jj & ~jj;                               % Markov blanket
    ee    = ~bb & ~jj;                                % external states
    b     = find(bb);
    e     = find(ee);
    s     = b(find( any(L(b,e),2)));
    a     = b(find(~any(L(b,e),2)));
    
    % partition
    %----------------------------------------------------------------------
    x{1,i} = spm_cat(z(a));
    x{2,i} = spm_cat(z(s));
    x{3,i} = spm_cat(z(j));
    
    % states accounted for
    %----------------------------------------------------------------------
    nn   = nn | ~ee;
    
    % plot
    %----------------------------------------------------------------------
    if GRAPHICS
        plot3(X(a,1),X(a,2),X(a,3),'.r','MarkerSize',24), hold on
        plot3(X(s,1),X(s,2),X(s,3),'.m','MarkerSize',24), hold on
        plot3(X(j,1),X(j,2),X(j,3),'.b','MarkerSize',24), hold on
        axis image
    end
    
end






