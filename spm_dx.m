function [dx] = spm_dx(dfdx,f,t)
% returns dx(t) = (expm(dfdx*t) - I)*inv(dfdx)*f
% FORMAT [dx] = spm_dx(dfdx,f,[t])
% dfdx   = df/dx
% f      = dx/dt
% t      = integration time: (default t = Inf);
%          if t is a cell (i.e., {t}) then t is set to:
%          exp(t - log(diag(-dfdx))
%
% dx     = x(t) - x(0)
%--------------------------------------------------------------------------
% Integration of a dynamic system using local linearization.  This scheme
% accommodates nonlinearities in the state equation by using a functional of
% f(x) = dx/dt.  This uses the equality
%
%             expm([0   0     ]) = (expm(t*dfdx) - I)*inv(dfdx)*f
%                  [t*f t*dfdx]
%
% When t -> Inf this reduces to
%
%              dx(t) = -inv(dfdx)*f
%
% These are the solutions to the gradient ascent ODE
%
%            dx/dt   = k*f = k*dfdx*x =>
%
%            dx(t)   = expm(t*k*dfdx)*x(0)
%                    = expm(t*k*dfdx)*inv(dfdx)*f(0) -
%                      expm(0*k*dfdx)*inv(dfdx)*f(0)
%
% When f = dF/dx (and dfdx = dF/dxdx), dx represents the update from a
% Gauss-Newton ascent on F.  This can be regularised by specifying {t}
% A heavy regularization corresponds to t = -4 and a light
% regularization would be t = 4. This version of spm_dx uses an augmented
% system and the Pade approximation to compute requisite matrix
% exponentials
%
% references:
%
% Friston K, Mattout J, Trujillo-Barreto N, Ashburner J, Penny W. (2007).
% Variational free energy and the Laplace approximation. NeuroImage.
% 34(1):220-34
%
% Ozaki T (1992) A bridge between nonlinear time-series models and
% nonlinear stochastic dynamical systems: A local linearization approach.
% Statistica Sin. 2:113-135.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dx.m 7144 2017-07-31 13:55:55Z karl $

% defaults
%--------------------------------------------------------------------------
nmax  = 512;                        % threshold for numerical approximation
if nargin < 3, t = Inf; end         % integration time
xf    = f; f = spm_vec(f);          % vectorise
n     = length(f);                  % dimensionality

% t is a regulariser
%--------------------------------------------------------------------------
sw  = warning('off','MATLAB:log:logOfZero');
if iscell(t)
    
    % relative integration time
    %----------------------------------------------------------------------
    t      = t{:};
    if isscalar(t)
        t  = exp(t - spm_logdet(dfdx)/n);
    else
        t  = exp(t - log(diag(-dfdx)));
    end
    
end
warning(sw);

% use a [pseudo]inverse if all t > TOL
%==========================================================================
if min(t) > exp(16)
    
    dx = -spm_pinv(dfdx)*f;
    
else
    
    % ensure t is a scalar or matrix
    %----------------------------------------------------------------------
    if isvector(t), t = diag(t); end
    
    % augment Jacobian and take matrix exponential
    %======================================================================
    J = spm_cat({0  []      ;
                t*f t*dfdx});
    
    % solve using matrix expectation
    %----------------------------------------------------------------------
    if n  <= nmax
        dx    = spm_expm(J);
        dx    = dx(:,1);
    else       
        x     = sparse(1,1,1,n + 1,1);
        dx    = expv(1,J,x);
    end
    
    % recover update
    %----------------------------------------------------------------------
    dx    = dx(2:end);
    
end
dx = spm_unvec(real(dx),xf);

return


%==========================================================================
%  Roger B. Sidje (rbs@maths.uq.edu.au)
%  EXPOKIT: Software Package for Computing Matrix Exponentials.
%  ACM - Transactions On Mathematical Software, 24(1):130-156, 1998

function  [w, err, hump] = expv( t, A, v, tol, m )
%  FOTMAT [w, err, hump] = expv( t, A, v, tol, m )
%  EXPV computes an approximation of w = exp(t*A)*v for a
%  general matrix A using Krylov subspace  projection techniques.
%  It does not compute the matrix exponential in isolation but instead,
%  it computes directly the action of the exponential operator on the
%  operand vector. This way of doing so allows for addressing large
%  sparse problems. The matrix under consideration interacts only
%  via matrix-vector products (matrix-free method).
%
%  w = expv( t, A, v )
%  computes w = exp(t*A)*v using a default tol = 1.0e-7 and m = 30.
%
%  [w, err] = expv( t, A, v )
%  renders an estimate of the error on the approximation.
%
%  [w, err] = expv( t, A, v, tol )
%  overrides default tolerance.
%
%  [w, err, hump] = expv( t, A, v, tol, m )
%  overrides default tolerance and dimension of the Krylov subspace,
%  and renders an approximation of the `hump'.
%
%  The hump is defined as:
%          hump = max||exp(sA)||, s in [0,t]  (or s in [t,0] if t < 0).
%  It is used as a measure of the conditioning of the matrix exponential
%  problem. The matrix exponential is well-conditioned if hump = 1,
%  whereas it is poorly-conditioned if hump >> 1. However the solution
%  can still be relatively fairly accurate even when the hump is large
%  (the hump is an upper bound), especially when the hump and
%  ||w(t)||/||v|| are of the same order of magnitude (further details in
%  reference below).
%
%  Example 1:
%  ----------
%    n = 100;
%    A = rand(n);
%    v = eye(n,1);
%    w = expv(1,A,v);
%
%  Example 2:
%  ----------
%    % generate a random sparse matrix
%    n = 100;
%    A = rand(n);
%    for j = 1:n
%        for i = 1:n
%            if rand < 0.5, A(i,j) = 0; end;
%        end;
%    end;
%    v = eye(n,1);
%    A = sparse(A); % invaluable for a large and sparse matrix.
%
%    tic
%    [w,err] = expv(1,A,v);
%    toc
%
%    disp('w(1:10) ='); disp(w(1:10));
%    disp('err =');     disp(err);
%
%    tic
%    w_matlab = expm(full(A))*v;
%    toc
%
%    disp('w_matlab(1:10) ='); disp(w_matlab(1:10));
%    gap = norm(w-w_matlab)/norm(w_matlab);
%    disp('||w-w_matlab|| / ||w_matlab|| ='); disp(gap);
%
%  In the above example, n could have been set to a larger value,
%  but the computation of w_matlab will be too long (feel free to
%  discard this computation).
%
%  See also MEXPV, EXPOKIT.

%  Roger B. Sidje (rbs@maths.uq.edu.au)
%  EXPOKIT: Software Package for Computing Matrix Exponentials.
%  ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
%__________________________________________________________________________

[n,n] = size(A);
if nargin == 3,
    tol = 1.0e-7;
    m = min(n,30);
end;
if nargin == 4,
    m = min(n,30);
end;

anorm = norm(A,'inf');
mxrej = 10;  btol  = 1.0e-7;
gamma = 0.9; delta = 1.2;
mb    = m; t_out   = abs(t);
nstep = 0; t_new   = 0;
t_now = 0; s_error = 0;
rndoff= anorm*eps;

k1    = 2; xm = 1/m; normv = norm(v); beta = normv;
fact  = (((m+1)/exp(1))^(m+1))*sqrt(2*pi*(m+1));
t_new = (1/anorm)*((fact*tol)/(4*beta*anorm))^xm;
s     = 10^(floor(log10(t_new))-1); t_new = ceil(t_new/s)*s;
sgn   = sign(t); nstep = 0;

w     = v;
hump  = normv;
while t_now < t_out
    nstep = nstep + 1;
    t_step = min( t_out-t_now,t_new );
    V = zeros(n,m+1);
    H = zeros(m+2,m+2);
    
    V(:,1) = (1/beta)*w;
    for j = 1:m
        p = A*V(:,j);
        for i = 1:j
            H(i,j) = V(:,i)'*p;
            p = p-H(i,j)*V(:,i);
        end;
        s = norm(p);
        if s < btol,
            k1 = 0;
            mb = j;
            t_step = t_out-t_now;
            break;
        end;
        H(j+1,j) = s;
        V(:,j+1) = (1/s)*p;
    end;
    if k1 ~= 0,
        H(m+2,m+1) = 1;
        avnorm = norm(A*V(:,m+1));
    end;
    ireject = 0;
    while ireject <= mxrej,
        mx = mb + k1;
        F = expm(sgn*t_step*H(1:mx,1:mx));
        if k1 == 0,
            err_loc = btol;
            break;
        else
            phi1 = abs( beta*F(m+1,1) );
            phi2 = abs( beta*F(m+2,1) * avnorm );
            if phi1 > 10*phi2,
                err_loc = phi2;
                xm = 1/m;
            elseif phi1 > phi2,
                err_loc = (phi1*phi2)/(phi1-phi2);
                xm = 1/m;
            else
                err_loc = phi1;
                xm = 1/(m-1);
            end;
        end;
        if err_loc <= delta * t_step*tol,
            break;
        else
            t_step = gamma * t_step * (t_step*tol/err_loc)^xm;
            s = 10^(floor(log10(t_step))-1);
            t_step = ceil(t_step/s) * s;
            if ireject == mxrej,
                error('The requested tolerance is too high.');
            end;
            ireject = ireject + 1;
        end;
    end;
    mx = mb + max( 0,k1-1 );
    w = V(:,1:mx)*(beta*F(1:mx,1));
    beta = norm( w );
    hump = max(hump,beta);
    
    t_now = t_now + t_step;
    t_new = gamma * t_step * (t_step*tol/err_loc)^xm;
    s = 10^(floor(log10(t_new))-1);
    t_new = ceil(t_new/s) * s;
    
    err_loc = max(err_loc,rndoff);
    s_error = s_error + err_loc;
end;
err  = s_error;
hump = hump / normv;

return


function E = padm( A, p )
%  FORMAT E = padm( A, p )
%  PADM computes the matrix exponential exp(A) using the irreducible
%  (p,p)-degree rational Pade approximation to the exponential function.
%
%  E = padm( A )
%  p is internally set to 6 (recommended and generally satisfactory).
%
%  See also CHBV, EXPOKIT and the MATLAB supplied functions EXPM and EXPM1.

%  Roger B. Sidje (rbs@maths.uq.edu.au)
%  EXPOKIT: Software Package for Computing Matrix Exponentials.
%  ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
%__________________________________________________________________________

if nargin == 1, p = 6; end;
[n,n] = size(A);

% Pade coefficients (1-based instead of 0-based as in the literature)
%--------------------------------------------------------------------------
c(1) = 1;
for k = 1:p
    c(k+1) = c(k)*((p+1-k)/(k*(2*p+1-k)));
end;

% Scaling
%--------------------------------------------------------------------------
s = norm(A,'inf');
if s > 0.5,
    s = max(0,fix(log(s)/log(2))+2);
    A = 2^(-s)*A;
end;

% Horner evaluation of the irreducible fraction (see ref. above)
%--------------------------------------------------------------------------
I = eye(n);
A2 = A*A;
Q = c(p+1)*I;
P = c(p)*I;
odd = 1;
for k = p-1:-1:1,
    if odd == 1,
        Q = Q*A2 + c(k)*I;
    else
        P = P*A2 + c(k)*I;
    end;
    odd = 1-odd;
end;
if odd == 1
    Q = Q*A;
    Q = Q - P;
    E = -(I + 2*(Q\P));
else
    P = P*A;
    Q = Q - P;
    E = I + 2*(Q\P);
end;

% Squaring
%--------------------------------------------------------------------------
for k = 1:s,
    E = E*E;
end;

return





