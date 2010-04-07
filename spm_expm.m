function [x] = spm_expm(J,x)
% approximate matrix exponential using a Taylor expansion
% FORMAT [y] = spm_expm(J,x)
% FORMAT [y] = spm_expm(J)
% y          = expm(J)*x:
% y          = expm(J);
%
% This routine covers and extends expm  functionality  by  using  a
% comoutationally  expedient  approximation  that can handle sparse
% matrices when dealing with the special case of expm(J)*x, where x
% is a vector, in an efficient fashion
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_expm.m 3813 2010-04-07 19:21:49Z karl $



if nargin == 1 || nargin == 2

    % expm(J) use Pade approximation
    %======================================================================

    % ensure norm is < 1/2 by scaling by power of 2
    %----------------------------------------------------------------------
    J     = full(J);
    I     = eye(size(J));
    [f,e] = log2(norm(J,'inf'));
    s     = max(0,e + 1);
    J     = J/2^s;
    X     = J;
    c     = 1/2;
    E     = I + c*J;
    D     = I - c*J;
    q     = 6;
    p     = 1;
    for k = 2:q
        c   = c*(q - k + 1)/(k*(2*q - k + 1));
        X   = J*X;
        cX  = c*X;
        E   = E + cX;
        if p
            D = D + cX;
        else
            D = D - cX;
        end
        p = ~p;
    end
    
    E = D\E;  % E = inv(D)*E

    % Undo scaling by repeated squaring E = E^(2^s)
    %----------------------------------------------------------------------
    for k = 1:s
        E = E*E;
    end

    % Multiply by x if necessary
    %----------------------------------------------------------------------
    if nargin > 1
        x = E*x;
    else
        x = E;
    end
    return


% Alternative evaluation of exp(J)*x for large matrices
%==========================================================================
else

    % compute y = expm(J)*x = (1 + J + J*J/2! + J*J*J/3!  + ...)*x
    %----------------------------------------------------------------------
    fx    = J*x;
    s     = 0;
    for i = 2:64

        % accumulate high-order terms until convergence
        %-----------------------------------------------------------
        nx = norm(fx,1);
        if ~nx
            return
        end
        fx = fx/nx;

        % take care to accumulate scaling in log-space
        %------------------------------------------------------------------
        s(i,1) = s(i - 1) + log(nx) - log(i - 1);
        x(:,i) = fx;
        fx     = J*fx;

        % add expansion terms if convergence
        %------------------------------------------------------------------
        if s(i) < -16
            x   = x*exp(s);
            return
        end
    end

    % If no convergence
    %----------------------------------------------------------------------
    disp('reverting to full pade')
    x   = spm_expm(J)*x(:,1);

end
