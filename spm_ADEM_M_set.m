function [M] = spm_ADEM_M_set(M)
% sets indices and performs checks on hierarchical action models
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
%   M(i).pC = prior covariances of p model-parameters
%   M(i).hE = prior expectation of h hyper-parameters (input noise)
%   M(i).hC = prior covariances of h hyper-parameters (input noise)
%   M(i).gE = prior expectation of g hyper-parameters (state noise)
%   M(i).gC = prior covariances of g hyper-parameters (state noise)
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
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_ADEM_M_set.m 3333 2009-08-25 16:12:44Z karl $

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


% and priors covariances - p
%--------------------------------------------------------------------------
try
    M.pC;
catch
    % Assume fixed parameters
    %----------------------------------------------------------------------
    for i = 1:g
        p       = length(spm_vec(M(i).pE));
        M(i).pC = sparse(p,p);
    end
end

% check pC if user specified
%--------------------------------------------------------------------------
for i = 1:g
 
    % Assume fixed parameters if not specified
    %----------------------------------------------------------------------
    if isempty(M(i).pC)
        p       = length(spm_vec(M(i).pE));
        M(i).pC = sparse(p,p);
    end
 
    % convert variances to covariances if necessary
    %----------------------------------------------------------------------
    if length(M(i).pC) == 1
        M(i).pC = sparse(diag(M(i).pC));
    end
 
    % check size
    %----------------------------------------------------------------------
    if length(M(i).pC) ~= length(spm_vec(M(i).pE))
        errordlg(sprintf('please check: M(%i).pC',i))
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

% ensure action is specified
%--------------------------------------------------------------------------
for i = 1:g
    try
        a  = M(i).a;
    catch
        a  = sparse(0,0);
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
        M(i).f  = fcnchk(M(i).f,'x','v','a','P');
    end
    try
        f       = feval(M(i).f,x,v,a,M(i).pE);
        if length(spm_vec(x)) ~= length(spm_vec(f))
            errordlg(sprintf('please check nargout: M(%i).f(x,v,a,P)',i));
        end
    catch
        errordlg(sprintf('evaluation failure: M(%i).f(x,v,a,P)',i))
    end
 
    % check g(x,v,P)
    %----------------------------------------------------------------------
    try
        M(i).g = fcnchk(M(i).g,'x','v','a','P');
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
        errordlg(sprintf('evaluation failure: M(%i).g(x,v,a,P)',i))
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


% Hyperparameters and components (causes: Q V and hidden states R, W)
%==========================================================================
try, M.Q;  catch, M(1).Q  = []; end
try, M.R;  catch, M(1).R  = []; end
try, M.V;  catch, M(1).V  = []; end
try, M.W;  catch, M(1).W  = []; end
try, M.hE; catch, M(1).hE = []; end
try, M.gE; catch, M(1).gE = []; end

% check hyperpriors hE - [log]hyper-parameters and components
%--------------------------------------------------------------------------
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
    try, M(i).hC*M(i).hE; catch, M(i).hC = speye(length(M(i).hE))*256; end
    try, M(i).gC*M(i).gE; catch, M(i).gC = speye(length(M(i).gE))*256; end
    
    if isempty(M(i).hC), M(i).hC = speye(length(M(i).hE))*256; end
    if isempty(M(i).gC), M(i).gC = speye(length(M(i).gE))*256; end
    
    % check Q and R (precision components)
    %======================================================================

    
    % check components and assume i.i.d if not specified
    %----------------------------------------------------------------------
    if length(M(i).Q) > length(M(i).hE)
        M(i).hE = sparse(length(M(i).Q),1);
    end
    if length(M(i).Q) < length(M(i).hE)
        M(i).Q  = {speye(M(i).l,M(i).l)};
        M(i).hE = M(i).hE(1);
    end
    if length(M(i).R) > length(M(i).gE)
        M(i).gE = sparse(length(M(i).R),1);
    end
    if length(M(i).R) < length(M(i).gE)
        M(i).R  = {speye(M(i).n,M(i).n)};
        M(i).gE = M(i).gE(1);
    end
    
    % check consistency and sizes (Q)
    %----------------------------------------------------------------------
    for j = 1:length(M(i).Q)
        if length(M(i).Q{j}) ~= M(i).l
            errordlg(sprintf('wrong size; M(%d).Q{%d}',i,j))
        end
    end
    
    % check consistency and sizes (R)
    %----------------------------------------------------------------------
    for j = 1:length(M(i).R)
        if length(M(i).R{j}) ~= M(i).n
            errordlg(sprintf('wrong size; M(%d).R{%d}',i,j))
        end
    end
    
    % check V and W (expansion point for precisions)
    %======================================================================

    % check V and assume unit precision if improperly specified
    %----------------------------------------------------------------------
    if length(M(i).V) ~= M(i).l
        try
            M(i).V = speye(M(i).l,M(i).l)*M(i).V(1);
        catch
            M(i).V = speye(M(i).l,M(i).l);
        end
    end
    
    % remove fixed components if hyperparameters exist
    %----------------------------------------------------------------------
    if ~isempty(M(i).hE)
        M(i).V = sparse(M(i).l,M(i).l);
    end
                
    % check W and assume unit precision if improperly specified
    %----------------------------------------------------------------------
    if length(M(i).W) ~= M(i).n
        try
            M(i).W = speye(M(i).n,M(i).n)*M(i).W(1);
        catch
            M(i).W = speye(M(i).n,M(i).n);
        end
    end
    
    % remove fixed components if hyperparameters exist
    %----------------------------------------------------------------------
    if ~isempty(M(i).gE)
        M(i).W = sparse(M(i).n,M(i).n);
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
try M(1).E.dt; catch M(1).E.dt = 1; end
 
% embedding orders
%--------------------------------------------------------------------------
try M(1).E.d;  catch, if nx, M(1).E.d = 2;  else M(1).E.d = 0;  end, end
try M(1).E.n;  catch, if nx, M(1).E.n = 6;  else M(1).E.n = 0;  end, end

M(1).E.d = min(M(1).E.d,M(1).E.n);
 
% number of iterations
%--------------------------------------------------------------------------
try M(1).E.nD; catch, if nx, M(1).E.nD = 1;  else M(1).E.nD = 8; end, end
try M(1).E.nE; catch,        M(1).E.nE = 1;  end
try M(1).E.nM; catch,        M(1).E.nM = 8;  end
try M(1).E.nN; catch,        M(1).E.nN = 16; end
 
% checks on smoothness hyperparameter
%==========================================================================
for i = 1:g
    
    try, M(i).sv; catch, M(i).sv = M(1).E.s;   end
    try, M(i).sw; catch, M(i).sw = M(1).E.s;   end
        
    if ~isscalar(M(i).sv), M(i).sv = M(1).E.s; end
    if ~isscalar(M(i).sw), M(i).sw = M(1).E.s; end
end


% checks on estimability
%==========================================================================
 
% check that there are informative priors on the states or the causes
%--------------------------------------------------------------------------
Q     = ~norm(M(end).V,1);
for i = 1:(g - 1)
    P = norm(M(i).pC,1) > exp(8);
    if P & Q
        warndlg('please use informative priors on causes or parameters')
    end
end


