function [x,u,y] = spm_Markov_blanket(J,z,m,mj)
% FORMAT [x,u,y] = spm_Markov_blanket(J,z,m,mj)
% Markovian partition
% J  - Jacobian
% z  - {1 x N}  partition of states (indices)
% m  - number of internal states [default: 3]
%
% x  - {3 x n} particular partition of state indices
%     x{1,j} - active states of j-th partition
%     x{2,j} - sensory states of j-th partition
%     x{3,j} - internal states of j-th partition
%
% u  - location of partitions in scaling or embedding space
%
% y  - {3 x n} particular partition of partition indices
%     y{1,j} - active states of j-th partition
%     y{2,j} - sensory states of j-th partition
%     y{3,j} - internal states of j-th partition
%
% Partition or Grouping (coarse-scaling) operator
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_Markov_blanket.m 7655 2019-08-25 20:10:20Z karl $

% preliminaries
%--------------------------------------------------------------------------
GRAPHICS  = 1;                      % Graphics switch
nz        = length(z);              % number of partitions
if nargin < 3
    m = 3;                          % maximum size of internal states
end
if nargin < 4
   mj = ones(nz,1);                 % eligible internal states
end
if isempty(mj)
   mj = ones(nz,1);                 % eligible internal states
end


% Adjacency matrix (over z)
%--------------------------------------------------------------------------
for i = 1:nz
    for j = 1:nz
        Lij    = J(z{i},z{j});
        if any(any(Lij))
            L(i,j) = abs(norm(full(Lij)) > 1/128);
        else
            L(i,j) = 0;
        end
    end
end
L     = double(L);

% get Markov blanket
%--------------------------------------------------------------------------
B     = L + L' + L'*L;
B     = B - diag(diag(B));

% scaling space (defined by graph Laplacian)
%--------------------------------------------------------------------------
G     = L + L';
G     = G - diag(diag(G));
G     = G - diag(sum(G));
G     = expm(G);

% get principal dimensions of scaling space (u)
%--------------------------------------------------------------------------
if GRAPHICS
    [u,v] = eig(G,'nobalance');
    v     = abs(diag(v));
    for i = 1:nz
        [p,h] = hist(real(u(:,i)),16);
        dh    = h(2) - h(1) + exp(-16);
        p     = p(:)/sum(p)/dh;
        v(i)  = log(v(i)) - p'*log(p + exp(-16))*dh;
    end
    [v,j] = sort(real(v),'descend');
    u     = real(u(:,j));
end


% recursive (particular) partition into internal, sensory and active states
%--------------------------------------------------------------------------
nn    = zeros(nz,1);
for i = 1:nz
    
    % internal states (defined by graph Laplacian)
    %----------------------------------------------------------------------
    jj = ~(B*nn) & ~nn & mj;
    if any(jj)
        
        % find densely coupled internal states (using the graph Laplacian)
        %------------------------------------------------------------------
        [g,j] = max(diag(G).*jj);
        if m > 1
            g      = G(:,j);
            g(j)   = 0;
            g(~jj) = 0;
            [g,k]  = sort(g,'descend');
            try
                j = [j; k(1:m - 1)];
            end
        end

        jj    = sparse(j,1,1,size(L,1),1) & jj;         % internal states
        bb    = B*jj & ~jj & ~nn;                       % Markov blanket
        ee    =  ~bb & ~jj & ~nn;                       % external states
        b     = find(bb);
        e     = find(ee);
        s     = b(find( any(L(b,e),2)));
        a     = b(find(~any(L(b,e),2)));
        
        % indices of individual states in the i-th particle
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
            a  = a(find(j,1));
            aa = sparse(a,1,1,size(L,1),1);
            ss = (L*aa | L'*aa) & ~aa & ~nn;
            a  = find(aa);
            s  = find(ss);
            j  = [];
            
            % indices of individual states in the i-th particle
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
            ss = sparse(s(1),1,1,nz,1);
            ss = ss | B*ss & ~nn;
            s  = find(ss);
            a  = [];
            j  = [];
            
            % indices of individual states in the i-th particle
            %--------------------------------------------------------------
            x{1,i} = [];
            x{2,i} = spm_cat(z(s));
            x{3,i} = [];
            
            % states accounted for (nn)
            %--------------------------------------------------------------
            nn   = nn | ss;
        end
    end
    
    % indices of partitions (i.e., n-states) in the i-th particle
    %----------------------------------------------------------------------
    y{1,i} = a;
    y{2,i} = s;
    y{3,i} = j;
    
    % plot
    %----------------------------------------------------------------------
    if all(nn) && numel(u) > 1
        
        % remove isolated (internal) states
        %--------------------------------------------------------------
        j     = [];
        for n = 1:size(x,2)
            if any(x{1,n}) || any(x{2,n})
                j = [j,n];
            end
        end
        x  = x(:,j);
        y  = y(:,j);
        
        if GRAPHICS,clf
            
            % colours for different particles
            %--------------------------------------------------------------
            nx            = size(x,2);
            [col,bol,msz] = spm_MB_col(nx);
            
            % plot partitions in embedding space (which particle)
            %--------------------------------------------------------------
            subplot(3,2,3)
            for k = 1:nx
                plot(u(y{1,k},1),u(y{1,k},2),'.','color',bol{k},'MarkerSize',msz), hold on
                plot(u(y{2,k},1),u(y{2,k},2),'.','color',bol{k},'MarkerSize',msz), hold on
                plot(u(y{3,k},1),u(y{3,k},2),'.','color',col{k},'MarkerSize',msz), hold on
            end
            axis square
            title(sprintf('Particles [%i n-states]',nz),'Fontsize',16)
            
            
            % plot particles in embedding space (which sort of state)
            %--------------------------------------------------------------
            subplot(3,2,4)
            for k = 1:nx
                plot(u(y{1,k},1),u(y{1,k},2),'.r','MarkerSize',msz), hold on
                plot(u(y{2,k},1),u(y{2,k},2),'.m','MarkerSize',msz), hold on
                plot(u(y{3,k},1),u(y{3,k},2),'.b','MarkerSize',msz), hold on
            end
            axis square
            title(sprintf('Markov partition [%i particles]',nx),'Fontsize',16)
            
            
            % plot particles in three embedding dimensions
            %--------------------------------------------------------------
            subplot(3,2,2)
            try
                for k = 1:nx
                    plot3(u(y{1,k},1),u(y{1,k},2),u(y{1,k},3),'.r','MarkerSize',msz), hold on
                    plot3(u(y{2,k},1),u(y{2,k},2),u(y{2,k},3),'.m','MarkerSize',msz), hold on
                    plot3(u(y{3,k},1),u(y{3,k},2),u(y{3,k},3),'.b','MarkerSize',msz), hold on
                end
            catch
                for k = 1:nx
                    plot(u(y{1,k},1),u(y{1,k},2),'.r','MarkerSize',msz), hold on
                    plot(u(y{2,k},1),u(y{2,k},2),'.m','MarkerSize',msz), hold on
                    plot(u(y{3,k},1),u(y{3,k},2),'.b','MarkerSize',msz), hold on
                end
            end
            
            axis square
            title('Embedding space','Fontsize',16)
            rotate3d(gca,'on')
            
            
            % Jacobian (ordered by partition and type)
            %--------------------------------------------------------------
            j = spm_vec(x');
            k = spm_vec(x );
            subplot(3,2,5),imagesc(-log(abs(J(k,k)) + exp(-4))),axis square
            subplot(3,2,6),imagesc(-log(abs(J(j,j)) + exp(-4))),axis square
            
            
            % Colors
            %--------------------------------------------------------------
            nj   = spm_length(x);
            msz  = fix(16 + 128/nj);
            j    = 1:nj;
            k    = spm_unvec(j,x')';
            j    = spm_unvec(j,x);
            subplot(3,2,5),hold on
            for q = 1:nx
                plot(j{1,q},ones(size(x{1,q})),'.','color',bol{q},   'MarkerSize',msz)
                plot(j{2,q},ones(size(x{2,q})),'.','color',bol{q},   'MarkerSize',msz)
                plot(j{3,q},ones(size(x{3,q})),'.','color',col{q},   'MarkerSize',msz)
                plot(j{1,q},zeros(size(x{1,q})) + nj,'.','color','r','MarkerSize',msz)
                plot(j{2,q},zeros(size(x{2,q})) + nj,'.','color','m','MarkerSize',msz)
                plot(j{3,q},zeros(size(x{3,q})) + nj,'.','color','b','MarkerSize',msz)
            end
            title(sprintf('Jacobian (by %i particles)',nx),'Fontsize',16)

            subplot(3,2,6),hold on
            for q = 1:nx
                plot(k{1,q},ones(size(x{1,q})),'.','color',bol{q},   'MarkerSize',msz)
                plot(k{2,q},ones(size(x{2,q})),'.','color',bol{q},   'MarkerSize',msz)
                plot(k{3,q},ones(size(x{3,q})),'.','color',col{q},   'MarkerSize',msz)
                plot(k{1,q},zeros(size(x{1,q})) + nj,'.','color','r','MarkerSize',msz)
                plot(k{2,q},zeros(size(x{2,q})) + nj,'.','color','m','MarkerSize',msz)
                plot(k{3,q},zeros(size(x{3,q})) + nj,'.','color','b','MarkerSize',msz)
            end
            title('Jacobian (by type)','Fontsize',16)
            
        end
        break
    end
    
end

return
