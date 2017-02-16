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
% $Id: FEP_MB_demo.m 7020 2017-02-16 13:15:35Z karl $


% default settings (GRAPHICS sets movies)
%--------------------------------------------------------------------------
% rng('default')
GRAPHICS = 1;

% create an adjacency matrix or Jacobian based upon a lattice
%==========================================================================

% within blanket coupling (intrinsic)
%--------------------------------------------------------------------------
n     = 2;
m     = 3;
Jii   = spm_cat({[], eye(n,n),  spm_speye(n,m);
                 zeros(n,n), [], [];
                 [], spm_speye(m,n), (randn(m,m)/8 - eye(m,m))});


% between blanket coupling (extrinsic)
%--------------------------------------------------------------------------
Jij   = spm_cat({[], zeros(n,n),  zeros(n,m);
                 randn(n,n),  [], [];
                 [], zeros(m,n), zeros(m,m)});


% an ensemble of blankets
%--------------------------------------------------------------------------
D     = 2;
N     = 8;                     % size of lattice
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
J     = spm_cat(A);


% Markov blanket partition
%--------------------------------------------------------------------------
x   = spm_Markov_blanket(J);
j   = spm_vec(x');
k   = spm_vec(x);

% dimension reduction (eliminating internal states)
%--------------------------------------------------------------------------
[Jr xr] = spm_A_reduce(J,x);
jr  = spm_vec(xr');
kr  = spm_vec(xr);

xx  = spm_Markov_blanket(Jr);
jj  = spm_vec(xx');
kk  = spm_vec(xx);


% parameters
%--------------------------------------------------------------------------
spm_figure('getwin','Figure 1'); clf

subplot(3,2,1),imagesc(abs(J(k,k))),axis square
subplot(3,2,2),imagesc(abs(J(j,j))),axis square

subplot(3,2,3),imagesc(abs(Jr(kr,kr))),axis square
subplot(3,2,4),imagesc(abs(Jr(jr,jr))),axis square

subplot(3,2,5),imagesc(abs(Jr(kk,kk))),axis square
subplot(3,2,6),imagesc(abs(Jr(jj,jj))),axis square

% Markov blanket - parents, children, and parents of children
%==========================================================================
function [J,y] = spm_A_reduce(J,x)

n     = 1;                        % number of generalised eigenvectors
nx    = size(x,2);                % number of partitions
for i = 1:nx
    Jsa    = full(J(x{2,i},x{1,i}));
    Jas    = full(J(x{1,i},x{2,i}));
    [e,~]  = eig(Jsa,Jas);
    v{i}   = e(:,1:n);
    u{i}   = pinv(v{i});
end
for i = 1:nx
    for j = 1:nx
        Jsa    = full(J(x{2,i},x{1,j}));
        Jas    = full(J(x{1,i},x{2,j}));
        A{i,j} = spm_cat({[], u{i}*Jas*v{j};
                          u{i}*Jsa*v{j}, []});
    end
    y{1,i} = (1:n) + (i - 1)*2*n;
    y{2,i} = (1:n) + (i - 1)*2*n + n;
end
J     = spm_cat(A);



function [x] = spm_Markov_blanket(J)

% Adjacency matrix
%--------------------------------------------------------------------------
m     = 3;
L     = sparse(double(any(J,3)));

% internal states (defined by graph Laplacian)
%----------------------------------------------------------------------
G     = L - diag(diag(L));
G     = G - diag(sum(G));
G     = expm(G);

% get Markov blanket and divide into sensory and active states
%----------------------------------------------------------------------
B     = double((L + L' + L*L'));
B     = B - diag(diag(B));

% recursive partition
%--------------------------------------------------------------------------
for i = 1:128
    
    % internal states (defined by graph Laplacian)
    %----------------------------------------------------------------------
    [g,j] = max(max(G));
    if g == 0
        break
    end
    g     = G(:,j);
    [g,j] = sort(g,'descend');   
    j     = j(1:m);                                   % internal cluster

    jj    = sparse(j,1,1,size(L,1),1);                % internal states
    bb    = B*jj & (1 - jj);                          % Markov blanket
    ee    = 1 - bb - jj;                              % external states
    b     = find(bb);
    e     = find(ee);
    s     = b(find( any(L(b,e),2)));
    a     = b(find(~any(L(b,e),2)));
    
    % partition
    %----------------------------------------------------------------------
    x{1,i} = a;
    x{2,i} = s;
    x{3,i} = j;
    
    % remove the Markov blanket and proceed
    %----------------------------------------------------------------------
    G(:,b) = 0;
    G(:,j) = 0;

end






