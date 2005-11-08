function [J] = spm_diff(varargin)
% matrix high-order differentials
% FORMAT [dfdx] = spm_diff(f,x,...,n)
% 
% f      - [inline] function f(x{1},...)
% x      - input argument[s]
% n      - arguments to differentiate w.r.t.
%
% dfdx          - df/dx{i}                     ; n =  i
% dfdx{p}...{q} - df/dx{i}dx{j}(q)...dx{k}(p)  ; n = [i j ... k]
%
% - a cunning recursive routine
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_diff.m 279 2005-11-08 19:11:28Z karl $



% create inline object
%--------------------------------------------------------------------------
f     = fcnchk(varargin{1});
x     = varargin(2:(end - 1));
n     = varargin{end};
m     = n(end);
xm    = spm_vec(x{m});
dx    = 1e-4;
J     = cell(1,length(xm));

% proceed to derivatives
%==========================================================================
if length(n) == 1

    % dfdx
    %----------------------------------------------------------------------
    f0    = feval(f,x{:});
    for i = 1:length(J)
        xi       = x;
        xmi      = xm;
        xmi(i)   = xmi(i) + dx;
        xi{n}    = spm_unvec(xmi,xi{n});
        fi       = feval(f,xi{:});
        J{i}     = spm_dfdx(fi,f0,dx);
    end

else

    % dfdxdxdx....
    %----------------------------------------------------------------------
    f0    = spm_diff(f,x{:},n(1:end - 1));
    for i = 1:length(J)
        xi       = x;
        xmi      = xm;
        xmi(i)   = xmi(i) + dx;
        xi{m}    = spm_unvec(xmi,xi{m});
        fi       = spm_diff(f,xi{:},n(1:end - 1));
        J{i}     = spm_dfdx(fi,f0,dx);
    end
    return

end

% return numeric array if for first order derivatives
%--------------------------------------------------------------------------
if length(n) == 1

    % vectorise f
    %----------------------------------------------------------------------
    f  = spm_vec(f0);
    
    % if there are no arguments to differentiate w.r.t. ...
    %----------------------------------------------------------------------
    if ~length(xm)
        J = sparse(length(f),0);
        return
    end
    
    % if there are no arguments to differentiate
    %----------------------------------------------------------------------
    if ~length(f)
        J = sparse(0,length(xm));
        return
    end

    % if f is a scalar 
    %----------------------------------------------------------------------
    if length(f) == 1
        J = spm_unvec(spm_vec(J),x{n});
        return
    end
    
    % if x{n} is a scalar
    %----------------------------------------------------------------------
    if length(xm) == 1
        J = spm_cat(J);
        return
    end
    
    % else f and xm are vectors return numeric array
    %----------------------------------------------------------------------
    for i = 1:length(J)
        J{i} = spm_vec(J{i});
    end
    J    = spm_cat(J);

end


function dfdx = spm_dfdx(f,f0,dx)
% cell subtraction
%--------------------------------------------------------------------------
if iscell(f)
   dfdx  = f;
   for i = 1:length(f(:))
       dfdx{i} = spm_dfdx(f{i},f0{i},dx);
   end
else
    dfdx  = (f - f0)/dx;
end
