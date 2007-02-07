function [varargout] = spm_diff(varargin)
% matrix high-order numerical differentiiation
% FORMAT [dfdx] = spm_diff(f,x,...,n,[V])
%
% f      - [inline] function f(x{1},...)
% x      - input argument[s]
% n      - arguments to differentiate w.r.t.
%
% dfdx          - df/dx{i}                     ; n =  i
% dfdx{p}...{q} - df/dx{i}dx{j}(q)...dx{k}(p)  ; n = [i j ... k]
%
% V      - cell array of matices that allow for differentiation w.r.t.
% to a linear trasnformation of the parameters: i.e., returns
% 
% df/dy{i};    x = V{i}y{i};    V = dx(i)/dy(i)
%
% - a cunning recursive routine
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_diff.m 731 2007-02-07 14:31:41Z karl $

% create inline object
%--------------------------------------------------------------------------
f     = fcnchk(varargin{1});

% parse input arguments
%--------------------------------------------------------------------------
if iscell(varargin{end})
    x = varargin(2:(end - 2));
    n = varargin{end - 1};
    V = varargin{end};
elseif isnumeric(varargin{end})
    x = varargin(2:(end - 1));
    n = varargin{end};
    V = cell(1,length(x));
else
    error('improper call')
end

% check trasnform matrices V = dx.dy
%--------------------------------------------------------------------------
for i = 1:length(x)
    try
        V{i};
    catch
        V{i} = [];
    end
    if ~length(V{i}) && any(n == i);
        V{i} = speye(length(spm_vec(x{i})));
    end
end

% initialise
%--------------------------------------------------------------------------
m     = n(end);
xm    = spm_vec(x{m});
dx    = exp(-8);
J     = cell(1,size(V{m},2));

% proceed to derivatives
%==========================================================================
if length(n) == 1

    % dfdx
    %----------------------------------------------------------------------
    f0    = feval(f,x{:});
    for i = 1:length(J)
        xi    = x;
        xmi   = xm + V{m}(:,i)*dx;
        xi{m} = spm_unvec(xmi,x{m});
        fi    = feval(f,xi{:});
        J{i}  = spm_dfdx(fi,f0,dx);
    end
    
    % return numeric array for first order derivatives
    %======================================================================

    % vectorise f
    %----------------------------------------------------------------------
    f  = spm_vec(f0);

    % if there are no arguments to differentiate w.r.t. ...
    %----------------------------------------------------------------------
    if ~length(xm)
        J = sparse(length(f),0);

    % or there are no arguments to differentiate
    %----------------------------------------------------------------------
    elseif ~length(f)
        J = sparse(0,length(xm));
        
    % or output is not a column vector
    %----------------------------------------------------------------------
    elseif size(f0,2) > 1
        for i = 1:size(J,2)
            J{i} = spm_vec(J{i});
        end
    end

    % concatenate into a matrix
    %----------------------------------------------------------------------
    J = spm_cat(J);
    
    % un-vectorise if f is a scalar
    %----------------------------------------------------------------------
    if length(f) == 1
        J = spm_unvec(J(:),x{m});
    end
    
    % un-vectorise x if is a scalar
    %----------------------------------------------------------------------
    if length(xm) == 1
        J = spm_unvec(J(:),f0);
    end
    
    % assign ouput argument and return
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
else
    dfdx  = (f - f0)/dx;
end
