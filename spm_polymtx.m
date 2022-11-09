function [b,D,H,o] = spm_polymtx(x,K,FUN)
% Create basis functions for polynomial expansion
% FORMAT [b,D,H,o] = spm_polymtx(x,K,FUN)
%
% x{i}   - domain of expansion (sample points): i = 1,...,N
% b      - expansion b = [... x{i}.^p.*x{j}.^q ...]: p,q = 0,...,(K - 1)
% D{i}   - first derivatives for each dimension: dx{i}/db
% H{i,j} - second derivatives : dx{j}dx{i}/dbdb
% o      - vector of expansion orders
%__________________________________________________________________________
%
% spm_polymtx creates a matrix for a polynomial expansion of order K - 1.
% With a second output argument, spm_polymtx produces the derivatives.
%
% b is a large prod(numel(x{i}) x K^N matrix corresponding to the Kroneckor
% tensor product of each N-dimensional domain. This is useful for dealing
% with vectorised N-arrays.
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 1996-2022 Wellcome Centre for Human Neuroimaging


% defaults
%--------------------------------------------------------------------------
if nargin < 2, K   = 3; end
if nargin < 3, FUN = 'POLY'; end

% Kroneckor form for n-dimensional DCTs
%--------------------------------------------------------------------------
n  = numel(x);
if n > 1
    
    % Kroneckor form of basis b
    %----------------------------------------------------------------------
    b     = 1;
    B     = cell(n,1);
    for i = 1:n
        B{i} = spm_polymtx(x(i),K,FUN);
        b    = kron(B{i},b);
    end

    % order of expansion o
    %----------------------------------------------------------------------
    o     = (1:K) - 1;
    for i = 2:n
        o = repmat(o,1,K);
        o = [o; kron((1:K) - 1,ones(1,K^(i - 1)))];
    end
    
    % derivatives (Kroneckor form)
    %----------------------------------------------------------------------
    if nargout > 1
        d     = cell(n,1);
        D     = cell(n,1);
        for i = 1:n
            [~,c] = spm_polymtx(x(i),K,FUN);
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
                    [~,~,c] = spm_polymtx(x(u),K,FUN);
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
    k  = sum(o) < K;
    o  = o(:,k);
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
N     = 32;
o     = (1:K) - 1;
x     = full(x{1}(:));
b     = zeros(size(x,1),K,'like',x);
switch FUN
    case {'POLY'}
        for k = 1:K
            b(:,k) = x.^(k - 1)/factorial(k - 1);
        end
    case {'DCT'}
        for k = 1:K
            b(:,k) = cos(pi*x*(k - 1)/N);
        end
    otherwise
        disp('Unknown expansion')
end


% derivatives
%--------------------------------------------------------------------------
if nargout > 1
    D = zeros(size(x,1),K,'like',x);
    switch FUN
        case {'POLY'}
            for k = 2:K
                D(:,k) = (x.^(k - 2)*(k - 1))/factorial(k - 1);
            end
        case {'DCT'}
            for k = 1:K
                D(:,k) = -(pi*sin((pi*x*(k - 1))/N)*(k - 1))/N;
            end
        otherwise
            disp('Unknown expansion')
    end
end

% Hessian
%--------------------------------------------------------------------------
if nargout > 2
    H = zeros(size(x,1),K,'like',x);
    switch FUN
        case {'POLY'}
            for k = 3:K
                H(:,k) = (x.^(k - 3)*(k - 1)*(k - 2))/factorial(k - 1);
            end
        case {'DCT'}
            for k = 1:K
                H(:,k) = -(pi^2*cos((pi*x*(k - 1))/N)*(k - 1)^2)/N^2;
            end
        otherwise
            disp('Unknown expansion')
    end
end
