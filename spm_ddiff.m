function [varargout] = spm_ddiff(varargin)
% matrix high-order numerical differentiation (double stencil)
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
% q      - (char) flag to preclude default concatenation of dfdx
%
% dfdx          - df/dx{i}                     ; n =  i
% dfdx{p}...{q} - df/dx{i}dx{j}(q)...dx{k}(p)  ; n = [i j ... k]
%
%
% This routine has the same functionality as spm_diff, however it
% uses two sample points to provide more accurate numerical (finite)
% differences that accommodate nonlinearities:
%
% dfdx  = (4*f(x + dx) - f(x + 2*dx) - 3*f(x))/(2*dx)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ddiff.m 7143 2017-07-29 18:50:38Z karl $

% step size for numerical derivatives
%--------------------------------------------------------------------------
global GLOBAL_DX
if ~isempty(GLOBAL_DX)
    dx = GLOBAL_DX;
else
    dx = exp(-8);
end

% create inline object
%--------------------------------------------------------------------------
f     = spm_funcheck(varargin{1});

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
for i = 1:length(x)
    try
        V{i};
    catch
        V{i} = [];
    end
    if isempty(V{i}) && any(n == i);
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
        xi     = x;
        xii    = x;
        xi{m}  = spm_unvec(xm + V{m}(:,i)*dx,  x{m});
        xii{m} = spm_unvec(xm + V{m}(:,i)*2*dx,x{m});
        J{i}   = spm_dfdx(f(xi{:}),f(xii{:}),f0,dx);
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
    varargout{1} = J;
    varargout{2} = f0;
    
else
    
    % dfdxdxdx....
    %----------------------------------------------------------------------
    f0        = cell(1,length(n));
    [f0{:}]   = spm_diff(f,x{:},n(1:end - 1),V);
    p         = true;
    
    for i = 1:length(J)
        xi     = x;
        xii    = x;
        xmi    = xm + V{m}(:,i)*dx;
        xmii   = xm + V{m}(:,i)*2*dx;
        xi{m}  = spm_unvec(xmi, x{m});
        xii{m} = spm_unvec(xmii,x{m});
        fi     = spm_diff(f,xi{:}, n(1:end - 1),V);
        fii    = spm_diff(f,xii{:},n(1:end - 1),V);
        J{i}   = spm_dfdx(fi,fii,f0{1},dx);
        p      = p & isnumeric(J{i});
    end
    
    % or differentiation of a scalar or vector
    %----------------------------------------------------------------------
    if p && q
        J = spm_dfdx_cat(J);
    end
    varargout = [{J} f0];
end


function dfdx = spm_dfdx(f,ff,f0,dx)
% cell subtraction
%__________________________________________________________________________
if iscell(f)
    dfdx  = f;
    for i = 1:length(f(:))
        dfdx{i} = spm_dfdx(f{i},ff{i},f0{i},dx);
    end
elseif isstruct(f)
    dfdx  = (4*spm_vec(f) - spm_vec(ff) - 3*spm_vec(f0))/(2*dx);
else
    dfdx  = (4*f - ff - 3*f0)/(2*dx);
end

return

function J = spm_dfdx_cat(J)
% concatenate into a matrix
%--------------------------------------------------------------------------
if isvector(J{1})
    if size(J{1},2) == 1
        J = spm_cat(J);
    else
        J = spm_cat(J')';
    end
end


