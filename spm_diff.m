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
% %W% Karl Friston %E%


% create inline object
%---------------------------------------------------------------------------
f     = fcnchk(varargin{1});
x     = varargin(2:(end - 1));
n     = varargin{end};
dx    = 1e-4;
J     = cell(1,length(x{n(end)}(:)));

% proceed to derivatives.
%---------------------------------------------------------------------------
if length(n) == 1

    % dfdx
    %-----------------------------------------------------------------------
    f0    = feval(f,x{:});
    for i = 1:length(x{n}(:))
        xi       = x;
        xi{n}(i) = xi{n}(i) + dx;
        fi       = feval(f,xi{:});
        J{i}     = spm_dfdx(fi,f0,dx);
    end

else

    % dfdxdxdx....
    %-----------------------------------------------------------------------
    f0    = spm_diff(f,x{:},n(1:end - 1));
    for i = 1:length(x{n(end)}(:))
        xi            = x;
        xi{n(end)}(i) = xi{n(end)}(i) + dx;
        fi            = spm_diff(f,xi{:},n(1:end - 1));
        J{i}          = spm_dfdx(fi,f0,dx);
    end
    return

end


try
    % if there are no arguments to differentiate w.r.t. ...
    %-----------------------------------------------------------------------
    if ~length(x{n}(:))
        J = sparse(length(f0),0);
    end
end

try
    % if there are no arguments to differentiate
    %-----------------------------------------------------------------------
    if ~size(J{1},1)
        J = sparse(0,length(x{n}(:)));
    end
end

try
    % reformat as a matrix
    %-----------------------------------------------------------------------
    if size(J{1},2) == 1
        J = spm_cat(J);
    end
end


function dfdx = spm_dfdx(f,f0,dx)
% cell subtraction
%--------------------------------------------------------------------------
if iscell(f)
   for i = 1:length(f)
       dfdx{i} = spm_dfdx(f{i},f0{i},dx);
   end
else
    dfdx  = (f - f0)/dx;
end
