function [f] = spm_DEM_basis(x,v,P)
% evaluates a parameterized set of basis functions
% problem
% FORMAT [f] = spm_DEM_basis(x,v,P)
%
% x   - hidden states
% v   - causal inputs
% P   - parameters
%
% returns f = sum(P(i)*B(x,i))
%
% where B(x,i) are basis functions
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_DEM_basis.m 3122 2009-05-13 16:05:59Z karl $
 
% basis set
%--------------------------------------------------------------------------
try, basis; catch, basis = 'radial'; end
 
% evaluate basis functions
%==========================================================================
switch basis
    
    case{'radial'}
        
        % Gaussian basis functions
        %------------------------------------------------------------------
        X = linspace(-2,2,length(P));
        W = 4*log(2)/(X(2) - X(1))^2;
        for i = 1:length(P)
            p(:,i) = exp(-W*(x - X(i)).^2);
        end
        f = (p*P(:))./sum(p,2);
        
    case{'spline'}
        
        % Natural spline
        %------------------------------------------------------------------
        X = linspace(-2,2,length(P));
        f = spline(X,P,x);
        
    case{'poly'}
        
        % Polynomial basis set
        %------------------------------------------------------------------
        f     = 0;
        for i = 1:length(P.p)
            B = x.^(i - 1);
            f = f + P.p(i)*B;
        end
        
    case{'dct'}
        
        % Discrete cosine basis set
        %------------------------------------------------------------------
        f     = 0;
        for i = 1:length(P.p)
            B = cos(i*pi*x/2);
            f = f + P.p(i)*B;
        end
end
 
f = spm_vec(f);



