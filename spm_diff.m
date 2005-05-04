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
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_diff.m 112 2005-05-04 18:20:52Z john $



% create inline object
%---------------------------------------------------------------------------
f     = fcnchk(varargin{1});
x     = varargin(2:(end - 1));
n     = varargin{end};
dx    = 1e-4;
J     = cell(size(x{n(end)}));

% proceed to derivatives.
%---------------------------------------------------------------------------
if length(n) == 1

    % dfdx
    %-----------------------------------------------------------------------
    f0    = feval(f,x{:});
    for i = 1:length(J(:))
        xi       = x;
        xi{n}(i) = xi{n}(i) + dx;
        fi       = feval(f,xi{:});
        J{i}     = spm_dfdx(fi,f0,dx);
    end

else

    % dfdxdxdx....
    %-----------------------------------------------------------------------
    f0    = spm_diff(f,x{:},n(1:end - 1));
    for i = 1:length(J(:))
        xi            = x;
        xi{n(end)}(i) = xi{n(end)}(i) + dx;
        fi            = spm_diff(f,xi{:},n(1:end - 1));
        J{i}          = spm_dfdx(fi,f0,dx);
    end
    return

end

% reformat as a matrix if appropriate 
%--------------------------------------------------------------------------
if length(n) == 1

    % if there are no arguments to differentiate w.r.t. ...
    %-----------------------------------------------------------------------
    if ~length(x{n})
        J = sparse(length(f0(:)),0);
    
    % if there are no arguments to differentiate
    %-----------------------------------------------------------------------
    elseif ~length(f0)
        J = sparse(0,length(x{n}(:)));
    end

    % if f and x{n} are vectors 
    %----------------------------------------------------------------------
    if min(size(f0)) == 1 & min(size(x{n})) == 1
        J = spm_cat(J(:)');

    % if f or x{n} are scalars 
    %----------------------------------------------------------------------
    elseif length(f0) == 1 | length(x{n}) == 1
        J = spm_cat(J);
    end
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
