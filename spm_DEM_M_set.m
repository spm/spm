function [M] = spm_DEM_M_set(M)
% sets indices and performs checks on hierarchical models
% FORMAT [M] = spm_DEM_M_set(M)
%
% for each level (i); required fields
%
%   M(i).g  = y(t)  = g(x,v,P)    {inline function, string or m-file}
%   M(i).f  = dx/dt = f(x,v,P)    {inline function, string or m-file}
%
% and
%
%   M(i).m  = number of inputs v(i + 1);
%   M(i).n  = number of states x(i);
%   M(i).l  = number of output v(i);
%
% or
%
%   M(i).x  = hidden states;
%   M(i).v  = causal states;
%
% for each level (i); optional fields
%
%   M(i).pE = prior expectation of p model-parameters
%   M(i).pC = prior covariances of p model-parameters
%   M(i).hE = prior expectation of h log-precision (cause noise)
%   M(i).hC = prior covariances of h log-precision (cause noise)
%   M(i).gE = prior expectation of g log-precision (state noise)
%   M(i).gC = prior covariances of g log-precision (state noise)
%   M(i).xC = prior covariances of states
%   M(i).Q  = precision components (input noise)
%   M(i).R  = precision components (state noise)
%   M(i).V  = fixed precision (input noise)
%   M(i).W  = fixed precision (state noise)
%
%
% sets fields, checks internal consistency of model specification and sets
% estimation parameters.  If a single hyperparameter is supplied i.i.d
% components are assumed (i.e., Q = I, R = I)
%--------------------------------------------------------------------------
%
%   M(1).E.s;     = smoothness (s.d. in time bins)
%   M(1).E.d;     = embedding order q(v)  (i.e., number of derivatives)
%   M(1).E.n;     = embedding order q(x)
%
% If the highest level involves any dynamic or static transformation
% of its inputs a further level is added with flat priors
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_DEM_M_set.m 5708 2013-10-22 09:20:59Z karl $

% order
%--------------------------------------------------------------------------
g      = length(M);
 
% set missing fields
%==========================================================================
 
% check for specification of hidden states
%--------------------------------------------------------------------------
if isfield(M,'f') && ~isfield(M,'n') && ~isfield(M,'x')
    msgbox('please specify hidden states or their number')
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
    [M.f] = deal(inline('sparse(0,1)','x','v','P'));
    [M.x] = deal(sparse(0,1));
    [M.n] = deal(0);
end
for i  = 1:g
    try
        fcnchk(M(i).f);
    catch
        M(i).f = inline('sparse(0,1)','x','v','P');
        M(i).x = sparse(0,1);
        M(i).n = 0;
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


% and priors covariances - p
%--------------------------------------------------------------------------
try
    M.pC;
catch
    
    % Assume fixed (zero variance) parameters
    %----------------------------------------------------------------------
    for i = 1:g
        p       = length(spm_vec(M(i).pE));
        M(i).pC = sparse(p,p);
    end
end

% check pC, if user specified
%--------------------------------------------------------------------------
for i = 1:g
    
    % number of parameters
    %----------------------------------------------------------------------
    np  = length(spm_vec(M(i).pE));
 
    % Assume fixed parameters if not specified
    %----------------------------------------------------------------------
    if isempty(M(i).pC)
        M(i).pC = sparse(np,np);
    end
 
    % convert variances to covariances if necessary
    %----------------------------------------------------------------------
    if isvector(M(i).pC)
        M(i).pC = sparse(diag(M(i).pC));
    end
    
    % convert variance to covariances if necessary
    %----------------------------------------------------------------------
    if isscalar(M(i).pC)
        M(i).pC = speye(np,np)*M(i).pC;
    end
 
    % check size
    %----------------------------------------------------------------------
    if length(M(i).pC) ~= np
        error('please check: M(%i).pC',i)
    end
 
end


% get inputs
%--------------------------------------------------------------------------
try
    v  = M(g).v;
catch
    v  = sparse(0,0);
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
        M(i).f  = fcnchk(M(i).f,'x','v','P');
    end
    try
        f       = feval(M(i).f,x,v,M(i).pE);
        if length(spm_vec(x)) ~= length(spm_vec(f))
            error('please check: M(%i).f(x,v,P)',i)
        end
 
    catch
        sprintf('??? evaluation failure: M(%i).f(x,v,P)',i)
        error(lasterror)
    end
    try M(i).fx = fcnchk(M(i).fx,'x','v','P'); end
    try M(i).fv = fcnchk(M(i).fv,'x','v','P'); end
    try M(i).fp = fcnchk(M(i).fp,'x','v','P'); end

 
    % check g(x,v,P)
    %----------------------------------------------------------------------
    try
        M(i).g = fcnchk(M(i).g,'x','v','P');
    end
    try
        M(i).m = length(spm_vec(v));
        v      = feval(M(i).g,x,v,M(i).pE);
        M(i).l = length(spm_vec(v));
        M(i).n = length(spm_vec(x));
 
        M(i).v = v;
        M(i).x = x;
 
    catch
        sprintf('??? evaluation failure: M(%i).g(x,v,P)',i)
        error(lasterror)
    end
    try M(i).gx = fcnchk(M(i).gx,'x','v','P'); end
    try M(i).gv = fcnchk(M(i).gv,'x','v','P'); end
    try M(i).gp = fcnchk(M(i).gp,'x','v','P'); end
    
end
    
% full priors on states
%--------------------------------------------------------------------------
try, M.xP; catch, M(1).xP = []; end
try, M.vP; catch, M(1).vP = []; end
for i = 1:g
    
    % hidden states
    %----------------------------------------------------------------------
    if isvector(M(i).xP), M(i).xP = diag(M(i).xP); end
    if length(M(i).xP) ~= M(i).n
        try
            M(i).xP = speye(M(i).n,M(i).n)*M(i).xP(1);
        catch
            M(i).xP = sparse(M(i).n,M(i).n);
        end
    end
    
    % hidden states
    %----------------------------------------------------------------------
    if isvector(M(i).vP), M(i).vP = diag(M(i).vP); end
    if length(M(i).vP) ~= M(i).l
        try
            M(i).vP = speye(M(i).l,M(i).l)*M(i).vP(1);
        catch
            M(i).vP = sparse(M(i).l,M(i).l);
        end
    end
end

% number of x (hidden states)
%--------------------------------------------------------------------------
nx     = sum(spm_vec(M.n));


% Hyperparameters and components (causes: Q V and hidden states R, W)
%==========================================================================
try, M.Q;  catch, M(1).Q  = []; end
try, M.R;  catch, M(1).R  = []; end
try, M.V;  catch, M(1).V  = []; end
try, M.W;  catch, M(1).W  = []; end
try, M.hE; catch, M(1).hE = []; end
try, M.gE; catch, M(1).gE = []; end
try, M.ph; catch, M(1).ph = []; end
try, M.pg; catch, M(1).pg = []; end

% check hyperpriors hE - [log]hyper-parameters and components
%--------------------------------------------------------------------------
pP    = 1;               % prior precision on log-precisions
for i = 1:g
    
    
    % make sure components are cell arrays
    %----------------------------------------------------------------------
    if ~isempty(M(i).Q) && ~iscell(M(i).Q), M(i).Q = {M(i).Q}; end
    if ~isempty(M(i).R) && ~iscell(M(i).R), M(i).R = {M(i).R}; end 
    
    % check hyperpriors
    %======================================================================
    
    % vectorise
    %----------------------------------------------------------------------
    M(i).hE = spm_vec(M(i).hE);
    M(i).gE = spm_vec(M(i).gE);
    
    % check hyperpriors (expectations)
    %----------------------------------------------------------------------
    if isempty(M(i).hE), M(i).hE = sparse(length(M(i).Q),1); end
    if isempty(M(i).gE), M(i).gE = sparse(length(M(i).R),1); end
    
    % check hyperpriors (covariances)
    %----------------------------------------------------------------------
    try, M(i).hC*M(i).hE; catch, M(i).hC = speye(length(M(i).hE))/pP; end
    try, M(i).gC*M(i).gE; catch, M(i).gC = speye(length(M(i).gE))/pP; end
    
    if isempty(M(i).hC), M(i).hC = speye(length(M(i).hE))/pP; end
    if isempty(M(i).gC), M(i).gC = speye(length(M(i).gE))/pP; end
    
    % check Q and R (precision components)
    %======================================================================

    
    % check components and assume i.i.d if not specified
    %----------------------------------------------------------------------
    if length(M(i).Q) > length(M(i).hE)
        M(i).hE = sparse(length(M(i).Q),1) + M(i).hE(1);
    end
    if length(M(i).Q) < length(M(i).hE)
        M(i).Q  = {speye(M(i).l,M(i).l)};
        M(i).hE = M(i).hE(1);
    end
    if length(M(i).hE) > length(M(i).hC)
        M(i).hC = speye(length(M(i).Q))*M(i).hC(1);
    end
    if length(M(i).R) > length(M(i).gE)
        M(i).gE = sparse(length(M(i).R),1) + M(i).gE(1);
    end
    if length(M(i).R) < length(M(i).gE)
        M(i).R  = {speye(M(i).n,M(i).n)};
        M(i).gE = M(i).gE(1);
    end
    if length(M(i).gE) > length(M(i).gC)
        M(i).gC = speye(length(M(i).R))*M(i).gC(1);
    end
    
    % check consistency and sizes (Q)
    %----------------------------------------------------------------------
    for j = 1:length(M(i).Q)
        if length(M(i).Q{j}) ~= M(i).l
            error('wrong size; M(%d).Q{%d}',i,j)
        end
    end
    
    % check consistency and sizes (R)
    %----------------------------------------------------------------------
    for j = 1:length(M(i).R)
        if length(M(i).R{j}) ~= M(i).n
            error('wrong size; M(%d).R{%d}',i,j)
        end
    end
    
    % check V and W (lower bound on precisions)
    %======================================================================

    % check V and assume unit precision if improperly specified
    %----------------------------------------------------------------------
    if isvector(M(i).V), M(i).V = diag(M(i).V); end
    if length(M(i).V) ~= M(i).l
        try
            M(i).V = speye(M(i).l,M(i).l)*M(i).V(1);
        catch
            if isempty(M(i).hE) && isempty(M(i).ph)
                M(i).V = speye(M(i).l,M(i).l);
            else
                M(i).V = sparse(M(i).l,M(i).l);
            end
        end
    end
    
                
    % check W and assume unit precision if improperly specified
    %----------------------------------------------------------------------
    if isvector(M(i).W), M(i).W = diag(M(i).W); end
    if length(M(i).W) ~= M(i).n
        try
            M(i).W = speye(M(i).n,M(i).n)*M(i).W(1);
        catch
            if isempty(M(i).gE) && isempty(M(i).pg)
                M(i).W = speye(M(i).n,M(i).n);
            else
                M(i).W = sparse(M(i).n,M(i).n);
            end
        end
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
try M(1).E.d;  catch, if nx, M(1).E.d = 2; else M(1).E.d = 0;  end, end
try M(1).E.n;  catch, if nx, M(1).E.n = 6; else M(1).E.n = 0;  end, end

M(1).E.d = min(M(1).E.d,M(1).E.n);
 
% number of iterations
%--------------------------------------------------------------------------
try M(1).E.nD; catch, if nx, M(1).E.nD = 1; else M(1).E.nD = 8; end, end
try M(1).E.nE; catch,        M(1).E.nE = 8; end
try M(1).E.nM; catch,        M(1).E.nM = 8; end
try M(1).E.nN; catch,        M(1).E.nN = 8; end
 
% checks on smoothness hyperparameter
%==========================================================================
try, M = rmfield(M,'sv'); end
try, M = rmfield(M,'sw'); end
for i = 1:g
    
    try, M(i).sv; catch, M(i).sv = M(1).E.s;   end
    try, M(i).sw; catch, M(i).sw = M(1).E.s;   end
        
    if ~isscalar(M(i).sv), M(i).sv = M(1).E.s; end
    if ~isscalar(M(i).sw), M(i).sw = M(1).E.s; end
end

% check on linear approximation scheme
%==========================================================================
try
    M(1).E.linear;
catch
    M(1).E.linear = 0;
end

% checks on estimability
%==========================================================================
 
% check that there are informative priors on the states or the causes
%--------------------------------------------------------------------------
Q     = ~norm(M(end).V,1);
for i = 1:(g - 1)
    P = norm(M(i).pC,1) > exp(8);
    if P && Q
        warndlg('please use informative priors on causes or parameters')
    end
end

