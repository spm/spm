function [x] = spm_funfun(varargin)
% Utility function to evaluate functionals
% FORMAT [F] = spm_funfun({f1,x11,x12,..f2,x22,...)
%
% F     = f ... f2(f1(x11,x12,...),x22,...)) ... )
%
% e.g. spm_funfun(@(x) cos(x),2.1,@(x,a) x^a,2)
%
% which is cos(2.1)^2
%__________________________________________________________________________
% Copyright (C) 2005-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_funfun.m 5219 2013-01-29 17:07:07Z spm $

% iterate over functions
%--------------------------------------------------------------------------
while ~isempty(varargin)
    f = varargin{1};
    n = nargin(f);    
    try
        x        = {x,varargin{2:n}};
        varargin = varargin(n + 1:end);
    catch
        x        = varargin(2:(n + 1));
        varargin = varargin(n + 2:end);
    end
    x = feval(f,x{:});
end
