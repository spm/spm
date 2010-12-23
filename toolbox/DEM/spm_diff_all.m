function [varargout] = spm_diff_all(varargin)
% matrix high-order numerical differentiation
% FORMAT [dfdx] = spm_diff(f,x,...,n)
% FORMAT [dfdx] = spm_diff(f,x,...,n,V)
% FORMAT [dfdx] = spm_diff(f,x,...,n,'q')
%
% f      - [inline] function f(x{1},...)
% x      - input argument[s]
% n      - arguments to differentiate w.r.t.
%
% V      - cell array of matrices that allow for differentiation w.r.t.
% to a linear transformation of the parameters: i.e., returns
% 
% df/dy{i};    x = V{i}y{i};    V = dx(i)/dy(i)
%
% q      - flag to preclude default concatenation of dfdx
%
% dfdx          - df/dx{i}                     ; n =  i
% dfdx{p}...{q} - df/dx{i}dx{j}(q)...dx{k}(p)  ; n = [i j ... k]
%
%
% - a cunning recursive routine
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_diff_all.m 4146 2010-12-23 21:01:39Z karl $
 
% create inline object
%--------------------------------------------------------------------------
f     = varargin{1};
 
% parse input arguments
%--------------------------------------------------------------------------
if iscell(varargin{end})
    x = varargin(2:(end - 2));
    n = varargin{end - 1};
    V = varargin{end};
    q = 1;
elseif isnumeric(varargin{end})
    x = varargin(2:(end - 1));
    n = varargin{end};
    V = cell(1,length(x));
    q = 1;
elseif ischar(varargin{end})
    x = varargin(2:(end - 2));
    n = varargin{end - 1};
    V = cell(1,length(x));
    q = 0;
else
    error('improper call')
end
 
% check transform matrices V = dxdy
%--------------------------------------------------------------------------
for i = 1:length(x(:))
    try
        V{i};
    catch
        V{i} = [];
    end
    if isempty(V{i}) && any(n == i);
      nx = size(x{1},2);
      V{i} = speye(length(spm_vec(x{i}))/nx);
    end
end
 
% initialise
%--------------------------------------------------------------------------
m     = n(end);
xm    = x{m};
dx    = exp(-8);
J0    = cell(1,size(V{m},2));
J     = sparse(size(V{m},2)*size(V{m},2),nx);
% proceed to derivatives
%==========================================================================
if length(n) == 1
 
    % dfdx
    %----------------------------------------------------------------------
    f0    = feval(f,x{:});
    for i = 1:length(J0)
        xi    = x;
        xmi   = xm + V{m}(:,i*ones(1,nx))*dx;
        xi{m} = spm_unvec(xmi,x{m});
        fi    = feval(f,xi{:});
        J0{i}  = spm_dfdx(fi,f0,dx);
    end

    J = sparse(cell2mat(J0'));
    % return numeric array for first-order derivatives
    %======================================================================
 
    % vectorise f
    %----------------------------------------------------------------------
    f  = spm_vec(f0);
   
    
    % assign output argument and return
    %----------------------------------------------------------------------
    varargout{1} = J;
    varargout{2} = f0;
 
else
 
    % dfdxdxdx....
    %----------------------------------------------------------------------
    f0        = cell(1,length(n));
    [f0{:}]   = spm_diff(f,x{:},n(1:end - 1),V);
    
    for i = 1:length(J)
        xi    = x;
        xmi   = xm + V{m}(:,i)*dx;
        xi{m} = spm_unvec(xmi,x{m});
        fi    = spm_diff(f,xi{:},n(1:end - 1),V);
        J{i}  = spm_dfdx(fi,f0{1},dx);
    end
    varargout = {J f0{:}};
end
 
 
function dfdx = spm_dfdx(f,f0,dx)
% cell subtraction
%--------------------------------------------------------------------------
if iscell(f)
    dfdx  = f;
    for i = 1:length(f(:))
        dfdx{i} = spm_dfdx(f{i},f0{i},dx);
    end
elseif isstruct(f)
    dfdx  = (spm_vec(f) - spm_vec(f0))/dx;
else
    dfdx  = (f - f0)/dx;
end
 
 
function is = isvec(v)
% isvector(v) returns true if v is 1-by-n or n-by-1 where n>=0
%__________________________________________________________________________
 
% vec if just two dimensions, and one (or both) unity
%--------------------------------------------------------------------------
is = length(size(v)) == 2 && isnumeric(v);
is = is && (size(v,1) == 1 || size(v,2) == 1);
