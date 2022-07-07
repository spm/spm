function [varargout] = spm_diff(varargin)
% Matrix high-order numerical differentiation
% FORMAT [dfdx] = spm_diff(f,x,...,n)
% FORMAT [dfdx] = spm_diff(f,x,...,n,V)
% FORMAT [dfdx] = spm_diff(f,x,...,n,'q')
%
% f      - [name or handle] function f(x{1},...)
% x      - input argument[s]
% n      - arguments to differentiate w.r.t.
%
% V      - cell array of matrices that allow for differentiation w.r.t.
%          to a linear transformation of the parameters: i.e., returns
%
%          df/dy{i};    x = V{i}y{i};    V = dx(i)/dy(i)
%
% q      - (char) flag to preclude default concatenation of dfdx
%
% dfdx          - df/dx{i}                     ; n =  i
% dfdx{p}...{q} - df/dx{i}dx{j}(q)...dx{k}(p)  ; n = [i j ... k]
%
%
% This routine has the same functionality as spm_ddiff, however it uses one
% sample point to approximate gradients with numerical (finite)
% differences:
%
% dfdx  = (f(x + dx)- f(x))/dx
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2003-2022 Wellcome Centre for Human Neuroimaging


% step size for numerical derivatives
%--------------------------------------------------------------------------
global GLOBAL_DX
if ~isempty(GLOBAL_DX)
    dx = GLOBAL_DX;
else
    dx = exp(-8);
end

% create function handle
%--------------------------------------------------------------------------
f     = spm_funcheck(varargin{1});

% parse input arguments
%--------------------------------------------------------------------------
if iscell(varargin{end})
    x = varargin(2:(end - 2));
    n = varargin{end - 1};
    V = varargin{end};
    q = true;
elseif isnumeric(varargin{end})
    x = varargin(2:(end - 1));
    n = varargin{end};
    V = cell(1,length(x));
    q = true;
elseif ischar(varargin{end})
    x = varargin(2:(end - 2));
    n = varargin{end - 1};
    V = cell(1,length(x));
    q = false;
else
    error('Improper call.')
end

% check transform matrices V = dxdy
%--------------------------------------------------------------------------
for i = 1:length(x)
    try
        V{i};
    catch
        V{i} = [];
    end
    if isempty(V{i}) && any(n == i)
        V{i} = speye(spm_length(x{i}));
    end
end

% initialise
%--------------------------------------------------------------------------
m     = n(end);
xm    = spm_vec(x{m});
J     = cell(1,size(V{m},2));

% proceed to derivatives
%==========================================================================
if length(n) == 1
    
    % dfdx
    %----------------------------------------------------------------------
    f0    = f(x{:});
    for i = 1:length(J)
        xi    = x;
        xi{m} = spm_unvec(xm + V{m}(:,i)*dx,x{m});
        J{i}  = spm_dfdx(f(xi{:}),f0,dx);
    end
    
    % return numeric array for first-order derivatives
    %======================================================================
    
    % vectorise f
    %----------------------------------------------------------------------
    f  = spm_vec(f0);
    
    % if there are no arguments to differentiate w.r.t. ...
    %----------------------------------------------------------------------
    if isempty(xm)
        J = sparse(length(f),0);
        
    % or there are no arguments to differentiate
    %----------------------------------------------------------------------
    elseif isempty(f)
        J = sparse(0,length(xm));
    end
    
    % differentiation of a scalar or vector
    %----------------------------------------------------------------------
    if isnumeric(f0) && iscell(J) && q
        J = spm_dfdx_cat(J);
    end
    
    % assign output argument and return
    %----------------------------------------------------------------------
    varargout = {J, f0};
    
else
    
    % dfdxdxdx....
    %----------------------------------------------------------------------
    f0        = cell(1,length(n));
    [f0{:}]   = spm_diff(f,x{:},n(1:end - 1),V);
    p         = true;
    
    for i = 1:length(J)
        xi    = x;
        xmi   = xm + V{m}(:,i)*dx;
        xi{m} = spm_unvec(xmi,x{m});
        fi    = spm_diff(f,xi{:},n(1:end - 1),V);
        J{i}  = spm_dfdx(fi,f0{1},dx);
        p     = p & isnumeric(J{i});
    end
    
    % or differentiation of a scalar or vector
    %----------------------------------------------------------------------
    if p && q
        J = spm_dfdx_cat(J);
    end
    
    % assign output argument and return
    %----------------------------------------------------------------------
    varargout = [{J}, f0];
    
end


%==========================================================================
% function dfdx = spm_dfdx(f,f0,dx)
%==========================================================================
function dfdx = spm_dfdx(f,f0,dx)
% numerical differences

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


%==========================================================================
% function J = spm_dfdx_cat(J)
%==========================================================================
function J = spm_dfdx_cat(J)
% concatenate into a matrix

if isvector(J{1})
    if size(J{1},2) == 1
        J = spm_cat(J);
    else
        J = spm_cat(J')';
    end
end
