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
%   M(i).pE = Prior expectation of p model parameters
%   M(i).pC = Prior covariances of p model parameters
%   M(i).hE = Prior expectation of h hyper-parameters
%   M(i).hC = Prior covariances of h hyper-parameters
%   M(i).Q  = free  precision components
%   M(i).V  = fixed precision
%
%
% sets fields, checks internal consistency of model specification and sets
% estimation parameters
%--------------------------
%
%   M(1).E.s;     = smoothness (seconds)
%   M(1).E.dt;    = time step
%   M(1).E.pt;    = prediction interval {L(t + pt)}
%   M(1).E.d;     = approximation order q(x,v)
%   M(1).E.n;     = order of embedding  (n >= d)
%
% If the highest level involves any dynamic or static transformation
% of its inputs a further level is added with flat priors
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id$

% order
%--------------------------------------------------------------------------
g      = length(M);
 
% set missing fields
%==========================================================================
 
% check for specifiaction of hidden states
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
    [M.f] = deal(inline('sparse(0,1)','x','v','P'));
    [M.x] = deal(sparse(0,1));
    [M.n] = deal(0);
end
 
% consistency and format check on states, parameters and functions
%==========================================================================


% prior expectation of parameters M.pE and initial values M.P
%--------------------------------------------------------------------------
try
    M.pE;
catch
    % Assume fixed parameters
    %----------------------------------------------------------------------
    for i = 1:(g - 1)
        M(i).pE = sparse(0,0);
    end
end
try
    M.P;
catch
    % Assume prior expectations
    %----------------------------------------------------------------------
    for i = 1:(g - 1)
        M(i).P = M(i).pE;
    end
end
for i = 1:(g - 1)
    if length(spm_vec(M(i).P)) ~= length(spm_vec(M(i).pE))
        errordlg(sprintf('please check: M(%i).pE/P',i))
    end
end

% and priors covariances - p
%--------------------------------------------------------------------------
try
    M.pC;
catch
    % Assume fixed parameters
    %----------------------------------------------------------------------
    for i = 1:(g - 1)
        p       = length(spm_vec(M(i).pE));
        M(i).pC = sparse(p,p);
    end
end

% check pC
%--------------------------------------------------------------------------
for i = 1:(g - 1)
 
    % Assume fixed parameters if not specified
    %----------------------------------------------------------------------
    if min(size(M(i).pC)) == 0
        p       = length(spm_vec(M(i).pE));
        M(i).pC = sparse(p,p);
    end
 
    % convert variances to covariances if necessary
    %----------------------------------------------------------------------
    if min(size(M(i).pC)) == 1
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
    v  = [];
end
if ~length(v)
    try
        v = sparse(M(g - 1).m,1);
    end
end
if ~length(v)
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
 
    % check f(x,v,P)
    %----------------------------------------------------------------------
    try
        M(i).f  = fcnchk(M(i).f,'x','v','P');
    end
    try
        f       = feval(M(i).f,x,v,M(i).pE);
        if ~length(spm_vec(x)) == length(spm_vec(f))
            str = sprintf('please check: M(%i).f(x,v,P)',i);
            msgbox(str)
            error(' ')
        end
 
    catch
        errordlg(sprintf('evaluation failure: M(%i).f(x,v,P)',i))
    end
 
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
        errordlg(sprintf('evaluation failure: M(%i).g(x,v,P)',i))
    end
end
 
% remove empty levels
%--------------------------------------------------------------------------
try
    g  = min(find(~cat(1,M.m)));
    M  = M(1:g);
catch
    errordlg('please specify number of variables')
end
 
 
% Precision components
%==========================================================================
 
% check Q{i} (free) and V (fixed)
%--------------------------------------------------------------------------
try
    M.Q;
catch
    M(1).Q = [];
end
try
    M.V;
catch
    M(1).V = [];
end
 
for i = 1:(g - 1)
 
    % make sure components are cell arrays
    %----------------------------------------------------------------------
    if length(M(i).Q) & ~iscell(M(i).Q)
        M(i).Q = {M(i).Q};
    end
 
    % make sure each level has a component
    %----------------------------------------------------------------------
    if ~length(M(i).V) & ~length(M(i).Q)
        M(i).Q = {speye(M(i).l,M(i).l)};
    end
 
    % check sizes
    %----------------------------------------------------------------------
    for j = 1:length(M(i).Q)
        if length(M(i).Q{j}) ~= M(i).l
            errordlg(sprintf('wrong size; M(%d).Q{%d}',i,j))
        end
    end
end
 
 
% check V (expansion point for covariances) and hyperparameters
%--------------------------------------------------------------------------
for i = 1:g
 
    % check V and assume V = 0 if improperly specified
    %----------------------------------------------------------------------
    if length(M(i).V) ~= M(i).l
        try
            M(i).V = speye(M(i).l,M(i).l)*M(i).V(1);
            warndlg(sprintf('expanding M(%i).V',i))
        catch
            M(i).V = sparse(M(i).l,M(i).l);
        end
    end
end
 
% Hyperparameters
%==========================================================================
 
% check hyperpriors h - [log]hyper-parameters and starting values M.h
%--------------------------------------------------------------------------
try
    M.hE;
catch
    for i = 1:g
        % Assume prior expectation of 0
        %----------------------------------------------------------------------
        M(i).hE = sparse(length(M(i).Q),1);
    end
end
try
    M.h;
catch
    % Assume prior expectations
    %----------------------------------------------------------------------
    for i = 1:g
        M(i).h = M(i).hE;
    end
end
for i = 1:g
 
    % check h and assume h = 0 if Q encodes a component (not in V)
    %----------------------------------------------------------------------
    M(i).hE = M(i).hE(:);
    M(i).h  = M(i).h(:);
    if length(M(i).hE) ~= length(M(i).Q)
        errordlg(sprintf('please check: M(%i).hE/Q',i))
    end
    if length(M(i).h) ~= length(M(i).hE)
        errordlg(sprintf('please check: M(%i).hE/h',i))
    end
end

% and prior covariances - h
%--------------------------------------------------------------------------
try
    M.hC;
catch
    % Assume fixed parameters
    %----------------------------------------------------------------------
    for i = 1:(g - 1)
        h       = length(M(i).Q);
        M(i).hC = speye(h,h)*16;
    end
end

% check size of hC
%--------------------------------------------------------------------------
for i = 1:(g - 1)
    if length(M(i).hC) ~= length(M(i).Q)
        errordlg(sprintf('please check: M(%i).hC',i))
    end
end

 
% estimation parameters M(1).E.s, n,...
%==========================================================================
% E.s;                               % smoothness (seconds)
% E.dt;                              % time step
% E.pt;                              % prediction interval {L(t + pt)}
% E.d;                               % approximation order of q(x,v)
% E.n;                               % order of embedding (n >= d)
nx     = sum(cat(1,M.n));            % number of x (hidden states)
 
% temporal smoothness - s.d. of kernel
%--------------------------------------------------------------------------
try M(1).E.s;  catch, if nx, M(1).E.s = 1; else M(1).E.s = 0;   end, end
 
% time step
%--------------------------------------------------------------------------
try M(1).E.dt; catch M(1).E.dt = 1;   end
try M(1).E.pt; catch M(1).E.pt = 1;   end
 
% embedding orders
%--------------------------------------------------------------------------
try M(1).E.d;  catch, if nx, M(1).E.d = 2;  else M(1).E.d = 1;  end, end
try M(1).E.n;  catch, if nx, M(1).E.n = 8;  else M(1).E.n = 1;  end, end
try M(1).E.r;  catch, M(1).E.r = 1;   end
 
% number of iterations
%--------------------------------------------------------------------------
try M(1).E.nD; catch, if nx, M(1).E.nD = 1; else M(1).E.nD = 8; end, end
try M(1).E.nE; catch, M(1).E.nE = 8;  end
try M(1).E.nM; catch, M(1).E.nM = 8;  end
try M(1).E.nI; catch, M(1).E.nI = 16; end
 
 
 
% checks on estimability
%==========================================================================
 
% check that there are informative priors on the states or the causes
%--------------------------------------------------------------------------
Q     = ~norm(M(end).V,1);
for i = 1:(g - 1)
    P = norm(M(i).pC,1) > 1024;
    if P & Q
        warndlg('please use informative priors on causes or parameters')
    end
end
