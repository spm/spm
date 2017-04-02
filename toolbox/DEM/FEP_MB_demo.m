function FEP_MB_demo
% This  routine illustrates a hierarchical decomposition of Markov blankets
% (of Markov blankets). It rests upon the dual operators of finding a
% partition (a Markov partition) and then using an adiabatic dimensional
% reduction (using the eigensolution of the Markov blanket). In brief,
% this means the states of particles at the next level become mixtures of
% the Markov blanket of particles at the level below.
%
% The ensuing hierarchical decomposition is illustrated in terms of
% Jacobian is and locations in a scathing space (evaluated using the graph
% Laplacian). This demonstration uses a fictive Jacobian that is created
% by hand.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: FEP_MB_demo.m 7051 2017-04-02 11:35:35Z karl $


% default settings
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
                 [], zeros(m,n),  zeros(m,m)});


% an ensemble of blankets
%--------------------------------------------------------------------------
D     = 2;                      % distance for seperation
N     = 8;                      % size of lattice
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
n     = [4 3 2 1];               % number of eigenvectors
m     = [3 3 3 3];               % number of internal states
for i = 1:N
    
    % Markov blanket partition
    %----------------------------------------------------------------------
    spm_figure('getwin',sprintf('Markov level %i',i));
    
    x{i} = spm_Markov_blanket(J{i},z{i},m(i));
    
    % dimension reduction (eliminating internal states)
    %----------------------------------------------------------------------
    if i < N
        [J{i + 1},z{i + 1}] = spm_A_reduce(J{i},x{i},n(i));
    end
    
end

return

% Markov blanket - parents, children, and parents of children
%==========================================================================
function [J,z] = spm_A_reduce(J,x,N)
% reduction of Markovian partition
% J  - Jacobian
% x  - {3 x N} indices of Markovian partition
% N  - number of eigenvectors to retain [default: 2]
%
% z  - {1 x N} indices of partition for the next level
%__________________________________________________________________________

% preliminaries
%--------------------------------------------------------------------------
nx    = size(x,2);                % number of partitions
if nargin < 3
    N = 2;                        % number of generalised eigenvectors
end

% reduction
%----------------------------------------------------------------------
for i = 1:nx
    Jii   = full(J(spm_vec(x(1:2,i)),spm_vec(x(1:2,i))));
    n     = min(N,size(Jii,1));
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



function [x,y] = spm_Markov_blanket(J,z,m)
% Markovian partition
% J  - Jacobian
% z  - {1 x N} indices of partition
% m  - number of internal states [default: 3]
%
% x  - {3 x N} indices of states of partitions
%     x{1,j} - active states of j-th partition
%     x{2,j} - sensory states of j-th partition
%     x{3,j} - internal states of j-th partition
%
% y  - {3 x N} indices of partition
%     y{1,j} - active states of j-th partition
%     y{2,j} - sensory states of j-th partition
%     y{3,j} - internal states of j-th partition
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
        L(i,j) = any(abs(Lij(:)) > 1e-16);
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
    [e,v] = eig(G,'nobalance');
    [v,j] = sort(real(diag(v)),'descend');
    try
        X = real(e(:,j(2:3)));
    catch
        X = real(e(:,j(1:2)));
    end
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
    j = ~(B*nn);
    if any(j)
        
        % find internal states
        %------------------------------------------------------------------
        [g,j] = max(max(G)'.*j);
        g     = G(:,j);
        [g,j] = sort(g,'descend');
        try
            j = j(1:m);                                 % internal states
        end
        
        jj    = sparse(j,1,1,size(L,1),1);              % internal states
        bb    = B*jj & ~jj;                             % Markov blanket
        ee    = ~bb & ~jj;                              % external states
        b     = find(bb);
        e     = find(ee);
        s     = b(find( any(L(b,e),2)));
        a     = b(find(~any(L(b,e),2)));
        
        % partition
        %------------------------------------------------------------------
        x{1,i} = spm_cat(z(a));
        x{2,i} = spm_cat(z(s));
        x{3,i} = spm_cat(z(j));
        
        % states accounted for (nn)
        %------------------------------------------------------------------
        nn   = nn | bb | jj;
        
    else
        
        % no internal states - find active states (not influenced by e)
        %------------------------------------------------------------------
        j = ~any(L(~nn,nn),2);
        if any(j)
            
            % sensory states connected with active states
            %--------------------------------------------------------------
            a  = find(~nn);
            a  = a(j(1));
            aa = sparse(a,1,1,size(L,1),1);
            ss = (L*aa | L'*aa) & ~aa & ~nn;
            a  = find(aa);
            s  = find(ss);
            j  = [];
            
            % partition
            %--------------------------------------------------------------
            x{1,i} = spm_cat(z(a));
            x{2,i} = spm_cat(z(s));
            x{3,i} = [];
            
            % states accounted for (nn)
            %--------------------------------------------------------------
            nn   = nn | aa | ss;
            
        elseif any(~nn)
            
            % sensory states connected with sensory states
            %--------------------------------------------------------------
            s  = find(~nn);
            ss = sparse(s(1),1,1,size(L,1),1);
            ss = (L*ss | L'*ss) & ~nn;
            s  = find(ss);
            a  = [];
            j  = [];
            
            % partition
            %--------------------------------------------------------------
            x{1,i} = [];
            x{2,i} = spm_cat(z(s));
            x{3,i} = [];
            
            % states accounted for (nn)
            %--------------------------------------------------------------
            nn   = nn | ss;
        end
    end
    
    % induces for the i-th particle
    %----------------------------------------------------------------------
    y{1,i} = a;
    y{2,i} = s;
    y{3,i} = j;
    
    % plot
    %----------------------------------------------------------------------
    if all(nn)
        if GRAPHICS,clf
            
            % plot different states
            %--------------------------------------------------------------
            subplot(3,2,1)
            nx    = size(x,2);
            for k = 1:nx
                plot(X(y{1,k},1),X(y{1,k},2),'.r','MarkerSize',24), hold on
                plot(X(y{2,k},1),X(y{2,k},2),'.m','MarkerSize',24), hold on
                plot(X(y{3,k},1),X(y{3,k},2),'.b','MarkerSize',24), hold on
            end
            axis image
            
            % plot different particles
            %--------------------------------------------------------------
            
            subplot(3,2,2)
            for k = 1:nx
                
                bol{k} = spm_softmax(log(rand(3,1))*2);
                col{k} = bol{k}*(1 - 1/2) + ones(3,1)/2;
                
                plot(X(y{1,k},1),X(y{1,k},2),'.','color',bol{k},'MarkerSize',24), hold on
                plot(X(y{2,k},1),X(y{2,k},2),'.','color',bol{k},'MarkerSize',24), hold on
                plot(X(y{3,k},1),X(y{3,k},2),'.','color',col{k},'MarkerSize',24), hold on
            end
            axis image
            
            % Jacobian
            %--------------------------------------------------------------
            j = spm_vec(x');
            k = spm_vec(x );
            subplot(3,2,3),imagesc(abs(J(k,k))),axis square, hold on
            subplot(3,2,4),imagesc(abs(J(j,j))),axis square, hold on
            
            % Colors
            %--------------------------------------------------------------
            nj   = spm_length(x);
            if nj > 32, msz  = 8; else, msz = 24; end
            j    = 1:nj;
            k    = spm_unvec(j,x')';
            j    = spm_unvec(j,x);
            subplot(3,2,3),hold on
            for q = 1:nx
                plot(j{1,q},ones(size(x{1,q})),'.','color',bol{q},   'MarkerSize',msz)
                plot(j{2,q},ones(size(x{2,q})),'.','color',bol{q},   'MarkerSize',msz)
                plot(j{3,q},ones(size(x{3,q})),'.','color',col{q},   'MarkerSize',msz)
                plot(j{1,q},zeros(size(x{1,q})) + nj,'.','color','r','MarkerSize',msz)
                plot(j{2,q},zeros(size(x{2,q})) + nj,'.','color','m','MarkerSize',msz)
                plot(j{3,q},zeros(size(x{3,q})) + nj,'.','color','b','MarkerSize',msz)
            end
            
            subplot(3,2,4),hold on
            for q = 1:nx
                plot(k{1,q},ones(size(x{1,q})),'.','color',bol{q},   'MarkerSize',msz)
                plot(k{2,q},ones(size(x{2,q})),'.','color',bol{q},   'MarkerSize',msz)
                plot(k{3,q},ones(size(x{3,q})),'.','color',col{q},   'MarkerSize',msz)
                plot(k{1,q},zeros(size(x{1,q})) + nj,'.','color','r','MarkerSize',msz)
                plot(k{2,q},zeros(size(x{2,q})) + nj,'.','color','m','MarkerSize',msz)
                plot(k{3,q},zeros(size(x{3,q})) + nj,'.','color','b','MarkerSize',msz)
            end   
            
        end
        
        break
    end
    
end






