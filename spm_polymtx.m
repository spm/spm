function [b,D,H,o] = spm_polymtx(x,K)
% Create basis functions for polynomial expansion
% FORMAT [b,D,H,o] = spm_polymtx(x,K)
%
% x{i}   - domain of expansion (sample points): i = 1,...,N
% b      - expansion b = [... x{i}.^p.*x{j}.^q ...]: p,q = 0,...,(K - 1)
% D{i}   - first derivatives for each dimension: dx{i}/db
% H{i,j} - second derivatives : dx{j}dx{i}/dbdb
% o      - vector of expansion orders
%__________________________________________________________________________
%
% spm_polymtx creates a matrix for a polynomial expansion of order K - 1.
% With a second output argument,spm_polymtx produces the derivatives.
%
% b is a large prod(numel(x{i}) x K^N matrix corresponding to the Kroneckor
% tensor product of each N-dimensional domain. This is useful for dealing
% with vectorised N-arrays.
%__________________________________________________________________________
% Copyright (C) 1996-2015 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_dctmtx.m 8000 2020-11-03 19:04:17Z karl $


% Kroneckor form for n-dimensional DCTs
%--------------------------------------------------------------------------
n  = numel(x);
if n > 1
    
    % Kroneckor form
    %----------------------------------------------------------------------
    b     = 1;
    o     = 0;
    O     = (1:K) - 1;
    B     = cell(n,1);
    for i = 1:n
        B{i} = spm_polymtx(x(i),K);
        b    = kron(B{i},b);
        o    = round(log(kron(exp(O),exp(o))));
    end
    
    % derivatives (Kroneckor form)
    %----------------------------------------------------------------------
    if nargout > 1
        d     = cell(n,1);
        D     = cell(n,1);
        for i = 1:n
            [~,c] = spm_polymtx(x(i),K);
            d{i}  = c;
            C     = 1;
            for u = 1:n
                if u == i
                    C = kron(d{u},C);
                else
                    C = kron(B{u},C);
                end
            end
            D{i}  = C;
        end
    end
    
    % Hessian (Kroneckor form)
    %----------------------------------------------------------------------
    if nargout > 2
        H     = cell(n,n);
        for i = 1:n
            
            % diagonal terms
            %--------------------------------------------------------------
            C     = 1;
            for u = 1:n
                if u == i
                    [~,~,c] = spm_polymtx(x(u),K);
                else
                    c = B{u};
                end
                C = kron(c,C);
            end
            H{i,i}  = C;
            
            % off-diagonal terms
            %--------------------------------------------------------------
            for j = (i + 1):n
                C     = 1;
                for u = 1:n
                    if u == i || u == j
                        C = kron(d{u},C);
                    else
                        C = kron(B{u},C);
                    end
                end
                H{i,j} = C;
                H{j,i} = C;
            end
        end
    end
    
    
    % remove high order terms
    %----------------------------------------------------------------------
    k  = o <= K;
    o  = o(k);
    b  = b(:,k);
    if nargout > 1
        for i = 1:n
            D{i} = D{i}(:,k);
        end
    end
    if nargout > 2
        for i = 1:n
            for j = 1:n
                H{j,i} = H{j,i}(:,k);
            end
        end
    end
    return
end

% polynomial expansion
%--------------------------------------------------------------------------
o     = (1:K) - 1;
x     = x{1}(:);
b     = zeros(size(x,1),K,'like',x);
for k = 1:K
    b(:,k) = x.^(k - 1);
end


% derivatives
%--------------------------------------------------------------------------
if nargout > 1
    D     = zeros(size(x,1),K,'like',x);
    for k = 2:K
        D(:,k) = (k - 1)*(x.^(k - 2));
    end
end

% Hessian
%--------------------------------------------------------------------------
if nargout > 2
    H     = zeros(size(x,1),K,'like',x);
    for k = 3:K
        H(:,k) = (k - 2)*(k - 1)*(x.^(k - 3));
    end
end



