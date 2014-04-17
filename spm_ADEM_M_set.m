function [M] = spm_ADEM_M_set(M)
% Set indices and perform checks on hierarchical action models
% FORMAT [M] = spm_ADEM_M_set(M)
%
% for each level (i); required fields
%
%   M(i).g  = y(t)  = g(x,v,a,P)    {inline function, string or m-file}
%   M(i).f  = dx/dt = f(x,v,a,P)    {inline function, string or m-file}
%
% and
%
%   M(i).m  = number of inputs v(i + 1);
%   M(i).n  = number of states x(i);
%   M(i).l  = number of output v(i);
%   M(i).k  = number of action a(i);
%
% or
%
%   M(i).x  = hidden states;
%   M(i).v  = causal states;
%   M(i).a  = action states;
%
% for each level (i); optional fields
%
%   M(i).pE = prior expectation of p model-parameters
%   M(i).V  = precision (input noise)
%   M(i).W  = precision (state noise)
%   M(i).U  = precision (action)
%
%
% sets fields, checks internal consistency of model specification and sets
% estimation parameters.  If (V,W) are not specified infinite precision is
% assumed.
%--------------------------------------------------------------------------
%
%   M(1).E.s;     = smoothness (s.d. in time bins)
%   M(1).E.d;     = embedding order q(v)  (i.e., number of derivatives)
%   M(1).E.n;     = embedding order q(x)
%
% If the highest level involves any dynamic or static transformation
% of its inputs a further level is added with flat priors
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_ADEM_M_set.m 5962 2014-04-17 12:47:43Z spm $


% order
%--------------------------------------------------------------------------
g      = length(M);
 
% set missing fields
%==========================================================================
 
% check for specification of hidden states
%--------------------------------------------------------------------------
if isfield(M,'f') && ~isfield(M,'n') && ~isfield(M,'x')
    msgbox('please specify hidden states or their number')
    error(' ')
end
 
 
% check supra-ordinate level and add one (with flat priors) if necessary
%--------------------------------------------------------------------------
try
    fcnchk(M(g).g);
    g      = g + 1;
    M(g).l = M(g - 1).m;
end
M(g).m = 0;
M(g).n = 0;
 
% default fields for static models (hidden states)
%--------------------------------------------------------------------------
if ~isfield(M,'f')
    [M.f] = deal(inline('sparse(0,1)','x','v','a','P'));
    [M.x] = deal(sparse(0,1));
    [M.n] = deal(0);
end
for i  = 1:g
    try
        fcnchk(M(i).f);
    catch
        M(i).f = inline('sparse(0,1)','x','v','a','P');
        M(i).x = sparse(0,1);
    end
end
 
% consistency and format check on states, parameters and functions
%==========================================================================
 
% prior expectation of parameters M.pE
%--------------------------------------------------------------------------
try
    M.pE;
catch
    % Assume fixed parameters
    %----------------------------------------------------------------------
    for i = 1:g
        M(i).pE = sparse(0,0);
    end
end

 
% get inputs
%--------------------------------------------------------------------------
try
    v  = M(g).v;
catch
    v  = sparse(0,1);
end
if isempty(v)
    try
        v = sparse(M(g - 1).m,1);
    end
end
if isempty(v)
    try
        v = sparse(M(g).l,1);
    end
end
M(g).l    = length(spm_vec(v));
M(g).v    = v;
 
% ensure action is specified
%--------------------------------------------------------------------------
for i = 1:g
    try
        a  = M(i).a;
    catch
        a  = sparse(0,1);
    end
    if isempty(a)
        try
            a = sparse(M(i).k,1);
        end
    end
    M(i).k    = length(spm_vec(a));
    M(i).a    = a;
end
 
 
% check functions
%--------------------------------------------------------------------------
for i = (g - 1):-1:1
    try
        x = M(i).x;
    catch
        x = sparse(M(i).n,1);
    end
    if isempty(x) && M(i).n
        x = sparse(M(i).n,1);
    end
 
    % check f(x,v,P)
    %----------------------------------------------------------------------
    try
        M(i).f = fcnchk(M(i).f);
        if nargin(M(i).f) ~= 4
            M(i).f = inline(char(M(i).f),'x','v','a','P');
        end
    end
    try
        f      = feval(M(i).f,x,v,a,M(i).pE);
        if length(spm_vec(x)) ~= length(spm_vec(f))
            errordlg(sprintf('please check nargout: G(%i).f(x,v,a,P)',i));
        end
    catch
        errordlg(sprintf('evaluation failure: G(%i).f(x,v,a,P)',i))
    end
 
    % check g(x,v,P)
    %----------------------------------------------------------------------
    try
        M(i).g = fcnchk(M(i).g);
        if nargin(M(i).g) ~= 4
            M(i).g = inline(char(M(i).g),'x','v','a','P');
        end
    end
    try
        M(i).m = length(spm_vec(v));
        v      = feval(M(i).g,x,v,a,M(i).pE);
        a      = M(i).a;
        M(i).k = length(spm_vec(a));
        M(i).l = length(spm_vec(v));
        M(i).n = length(spm_vec(x));
 
        M(i).a = a;
        M(i).v = v;
        M(i).x = x;
 
    catch
        errordlg(sprintf('evaluation failure: G(%i).g(x,v,a,P)',i))
    end
end
 
% remove empty levels
%--------------------------------------------------------------------------
try
    g  = find(~spm_vec(M.m),1);
    M  = M(1:g);
catch
    errordlg('please specify number of variables')
end
 
% number of x (hidden states)
%--------------------------------------------------------------------------
nx     = sum(spm_vec(M.n));
 
 
% check precisions
%==========================================================================
try, M.V;  catch, M(1).V  = []; end
try, M.W;  catch, M(1).W  = []; end

 
% check precisions
%--------------------------------------------------------------------------
for i = 1:g
 
    % check V and assume unit precision if improperly specified
    %----------------------------------------------------------------------
    if isvector(M(i).V), M(i).V = diag(M(i).V); end
    if length(M(i).V) ~= M(i).l
        try
            M(i).V = speye(M(i).l,M(i).l)*M(i).V(1);
        catch
            M(i).V = speye(M(i).l,M(i).l)*exp(32);
        end
    end
                
    % check W and assume unit precision if improperly specified
    %----------------------------------------------------------------------
    if isvector(M(i).W), M(i).W = diag(M(i).W); end
    if length(M(i).W) ~= M(i).n
        try
            M(i).W = speye(M(i).n,M(i).n)*M(i).W(1);
        catch
            M(i).W = speye(M(i).n,M(i).n)*exp(32);
        end
    end
end


% check restiction of precision for action (gain)
%==========================================================================
try, M(1).U;  catch, M(1).U = []; end

if isvector(M(1).U), M(1).U = diag(M(1).U); end
if length(M(1).U) ~= M(1).l
    try
        M(1).U = speye(M(1).l,M(1).l)*M(1).U(1);
    catch
        M(1).U = speye(M(1).l,M(1).l)*exp(2);
    end
end

 
% estimation parameters M(1).E.s, n,...
%==========================================================================
% E.s;                               % smoothness (seconds)
% E.dt;                              % time step
% E.d;                               % approximation order of q(x,v)
% E.n;                               % order of embedding (n >= d)
 
% temporal smoothness - s.d. of kernel
%--------------------------------------------------------------------------
try M(1).E.s;  catch, if nx, M(1).E.s = 1/2; else M(1).E.s = 0; end, end
 
% time step
%--------------------------------------------------------------------------
try M(1).E.dt; catch, M(1).E.dt = 1; end
 
% embedding orders
%--------------------------------------------------------------------------
try M(1).E.d;  catch, if nx, M(1).E.d = 2;  else M(1).E.d = 0;  end, end
try M(1).E.n;  catch, if nx, M(1).E.n = 6;  else M(1).E.n = 0;  end, end
 
M(1).E.d = min(M(1).E.d,M(1).E.n);

 
% checks on smoothness hyperparameter
%==========================================================================
for i = 1:g
    
    try, M(i).sv; catch,   M(i).sv = M(1).E.s; end
    try, M(i).sw; catch,   M(i).sw = M(1).E.s; end
        
    if ~isscalar(M(i).sv), M(i).sv = M(1).E.s; end
    if ~isscalar(M(i).sw), M(i).sw = M(1).E.s; end
end
