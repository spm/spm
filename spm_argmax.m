function [p f] = spm_argmax(fun,varargin)
% function minimisation using Gauss-Newton
% FORMAT [P f] = spm_argmax('fun',varargin,i)
%
% fun      - inline function f - fun(P,varargin)
% varargin - function arguments
% i        - argument to minimise: varargin{i}
%
% P   - optimised argument
% f   - optimised value of fun(P)
%
%--------------------------------------------------------------------------
% spm_argmax uses numerical derivatives and a and adaptive Gauss-Newton 
% scheme: see also spm_argmin.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_argmax.m 4297 2011-04-07 18:12:29Z karl $
 
 
% arguments
%--------------------------------------------------------------------------
x     = varargin(1:end - 1);
i     = varargin{end};

% Gauss-Newton scheme
%======================================================================
v     = -8;
for j = 1:64
    
    % time
    %----------------------------------------------------------------------
    tic
    
    % get derivatives of function w.r.t. argument
    %----------------------------------------------------------------------
    [dfdpp dfdp f] = spm_diff(fun,x{:},[i i]);
    dfdp           = spm_vec(dfdp);
    dfdpp          = spm_cat(dfdpp');
    dfdpp          = -spm_sqrtm(dfdpp'*dfdpp);
    
    % report change in f
    %----------------------------------------------------------------------
    try
        fprintf(' actual: %6.2e (%.2f sec)\n',full(f - P.f),toc)
    catch
        f0   = f;
        p    = spm_vec(x{i});
        P.f  = -Inf;
    end
    
    % If f is decreasing accept current estimates decrease regularization
    %----------------------------------------------------------------------
    if f > P.f
        P.x  = x;
        P.f  = f;
        v    = min(v + 1/8,4);
        v    = max(v,-8);
        str  = 'GN:(+)';
 
    else
        x    = P.x;
        v    = min(v - 2, -4);
        str  = 'GN:(-)';
    end

    % Gradient ascent
    %======================================================================
    dp   = spm_dx(dfdpp,dfdp,{v});
    df   = dfdp'*dp;
      
    % check stability
    %----------------------------------------------------------------------
    while df > 1e4
        v    = v - log(abs(df));
        dp   = spm_dx(dfdpp,dfdp,{v});
        df   = dfdp'*dp;
    end
    
    % update
    %---------------------------------------------------------------------- 
    p    = p + dp;
    x{i} = spm_unvec(p,x{i});
    
    % convergence
    %----------------------------------------------------------------------
    fprintf('%-6s: %i %6s %-6.2e %6s %-6.2e',str,j,'F:',full(P.f - f0),'dF predicted:',full(df))
    if j > 2 && df < 1e-6
        fprintf(' convergence\n')
        break
    end
end

% output
%----------------------------------------------------------------------
p = x{i};
