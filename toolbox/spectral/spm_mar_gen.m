function [v] = spm_mar_gen (w,A,C,n,ndisc)
% Generate data from MAR model
% FORMAT [v] = spm_mar_gen (w,A,C,n,ndisc)
%
% Generates n time steps of the MAR(p) process
%
%     v(k,:)' = w' + A1*v(k-1,:)' +...+ Ap*v(k-p,:)' + eta(k,:)', 
%
%  where A=[A1 ... Ap] is the coefficient matrix, and w is a vector of
%  intercept terms that is included to allow for a nonzero mean of the
%  process. The vectors eta(k,:) are independent Gaussian noise
%  vectors with mean zero and covariance matrix C.
%
%  This function is adapted from the ARFIT toolbox by Neumaier and
%  Schneider
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_mar_gen.m 1143 2008-02-07 19:33:33Z spm $

m       = size(C,1);                  % dimension of state vectors 
p       = size(A,2)/m;                % order of process

if (p ~= round(p)) 
    error('Bad arguments.'); 
end

if (length(w) ~= m | min(size(w)) ~= 1)
    error('Dimensions of arguments are mutually incompatible.')
end 
w       = w(:)';                      % force w to be row vector

% Discard the first ndisc time steps; if ndisc is not given as input
% argument, use default
if (nargin < 5) 
    ndisc = 10^3; 
end

% Compute Cholesky factor of covariance matrix C
R       = chol(C);                    % R is upper triangular

% Get ndisc+n independent Gaussian pseudo-random vectors with 
% covariance matrix C=R'*R
randvec = randn([ndisc+n,m])*R;

% Add intercept vector to random vectors
randvec = randvec + ones(ndisc+n,1)*w;

% Get transpose of system matrix A (use transpose in simulation because 
% we want to obtain the states as row vectors)
AT      = A';

% Take the p initial values of the simulation to equal the process mean, 
% which is calculated from the parameters A and w
if any(w)
    %  Process has nonzero mean    mval = inv(B)*w'    where 
    %             B = eye(m) - A1 -... - Ap; 
    %  Assemble B
    B    = eye(m);
    for j=1:p
        B = B - A(:, (j-1)*m+1:j*m);
    end
    %  Get mean value of process
    mval = w / B';
    
    %  The optimal forecast of the next state given the p previous
    %  states is stored in the vector x. The vector x is initialized
    %  with the process mean.
    x    = ones(p,1)*mval;
else
    %  Process has zero mean
    x    = zeros(p,m); 
end

% Initialize state vectors
u      = [x; zeros(ndisc+n,m)];

% Simulate n+ndisc observations. In order to make use of Matlab's
% vectorization capabilities, the cases p=1 and p>1 must be treated 
% separately.
if p==1
    for k=2:ndisc+n+1; 
        x(1,:) = u(k-1,:)*AT;
        u(k,:) = x + randvec(k-1,:);
    end
else
    for k=p+1:ndisc+n+p; 
        for j=1:p;
            x(j,:) = u(k-j,:)*AT((j-1)*m+1:j*m,:);
        end
        u(k,:) = sum(x)+randvec(k-p,:);
    end
end

% return only the last n simulated state vectors
v = u(ndisc+p+1:ndisc+n+p,:); 





