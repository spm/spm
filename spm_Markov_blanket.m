function [x,y] = spm_Markov_blanket(J,z,m,R)
% Markovian partition
% FORMAT [x,y] = spm_Markov_blanket(J,z,m,R)
% J  - Jacobian
% z  - {1 x N}  partition of states (indices)
% m  - number of internal states [default: 3]
%
% x  - {3 x n} particular partition of state indices
%     x{1,j} - active states of j-th partition
%     x{2,j} - sensory states of j-th partition
%     x{3,j} - internal states of j-th partition
%
% y  - {3 x n} particular partition of partition indices
%     y{1,j} - active states of j-th partition
%     y{2,j} - sensory states of j-th partition
%     y{3,j} - internal states of j-th partition
%
% Partition or Grouping (coarse-scaling) operator
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


% preliminaries
%--------------------------------------------------------------------------
if nargin < 3, m = 3;  end                % maximum size of internal states
if nargin < 4, R = []; end                % restiction matrix

% Adjacency matrix (over z)
%--------------------------------------------------------------------------
nz    = length(z);                        % number of partitions
for i = 1:nz
    for j = 1:nz
        Lij    = J(z{i},z{j});
        if any(any(Lij))
            L(i,j) = norm(full(Lij));
        else
            L(i,j) = 0;
        end
    end
end

% supress coupling (and apply restriction if specified)
%--------------------------------------------------------------------------
% if numel(R), L = L.*R; end
L(L < 1/64) = 0;

% get Markov blanket matrix
%--------------------------------------------------------------------------
B     = L + L' + L'*L;
B     = B - diag(diag(B));
B     = sparse(B);

% scaling space (defined by graph Laplacian)
%--------------------------------------------------------------------------
% G     = L + L';
% G     = G - diag(diag(G));
% G     = G - diag(sum(G));
% G     = expm(G);

% recursive (particular) partition into internal, sensory and active states
%--------------------------------------------------------------------------
nn    = zeros(nz,1);
for i = 1:nz
    
    % internal states (defined by graph Laplacian)
    %----------------------------------------------------------------------
    jj = ~(B*nn) & ~nn;
    ij = find(jj);
    if any(jj)
        
        % find densely coupled internal states (using the eigenmode of B)
        %------------------------------------------------------------------
        [v,s] = svds(B(ij,ij),1);
        [v,j] = sort(abs(v),'descend');
        j     = ij(j(1:min(m,numel(j))));
        
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
    
    % remove isolated (internal) states
    %----------------------------------------------------------------------
    if all(nn)
        j     = [];
        for n = 1:size(x,2)
            if any(x{1,n}) || any(x{2,n})
                j = [j,n];
            end
        end
        x  = x(:,j);
        y  = y(:,j);
        break
    end
    
end
