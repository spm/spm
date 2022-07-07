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

% Karl Friston
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging

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
