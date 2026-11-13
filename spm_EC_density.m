function [EC] = spm_EC_density(STAT,t,df,D)
% Returns the Euler characteristic (EC) density
% FORMAT function [EC] = spm_EC_density(STAT,t,df,[D])
%
% STAT  = 'Z','T','F',...
% t     = vector of statistical values
% df    = degrees of freedom [df1,df2]
% D     = dimension of SPM [default: 4]
%
% This routine returns the density of the expected Euler characteristic for
% a variety of commonly used statistics. It is based upon a generic
% expression for fields of the F statistic. This provides estimates for
% special cases: for example, with infinite degrees of freedom of the
% denominator, the F statistic behaves as a chi-squared statistic.
% Similarly, when the degrees of freedom of the numerator reduce to 1, the
% squared T statistic behaves as an F statistic. And as a Z statistic,
% when the denominator degrees of freedom tend to infinity.
%
% This routine returns densities over any arbitrary number of dimensions,
% D. The default is D = 4 dimensions returning a D + 1 densities, where the
% first corresponds to the zeroth dimension or point (i.e., using the usual
% cumulative density for the statistic in question).
%
% The current implementation is based upon an optimised functional form in
% Appendix B of Worsley et al (2004). The (polynomial) summations are
% grouped in a way to maximise efficiency. The gamma functions are
% evaluated in log space. This enables evaluation with very high degrees of
% freedom (with an enforced upper bound of 256).
%__________________________________________________________________________
%
% Reference : K.J. Worsley et al. / NeuroImage 23 (2004) S189–S195
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 1999-2022 Wellcome Centre for Human Neuroimaging

% set-up
%--------------------------------------------------------------------------
if nargin < 4, D = 4; end                         % dimension of SPM

% EC densities
%--------------------------------------------------------------------------
t   = t(:)';
EC  = zeros(D + 1, numel(t));
if STAT == 'Z'

    % Gaussian Field: SPM(Z)
    %----------------------------------------------------------------------
    for d = 1:(D + 1)
        if d > 1
            EC(d,:) = ECdF(d - 1,[1,Inf],t.^2).*(sign(t).^d)/2;
        else
            EC(d,:) = 1 - spm_Ncdf(t);
        end
    end

elseif STAT == 'T'

    % SPM(T)
    %----------------------------------------------------------------------
    for d = 1:(D + 1)
        if d > 1
            EC(d,:) = ECdF(d - 1,[1,df(2)],t.^2).*(sign(t).^d)/2;
        else
            EC(d,:) = 1 - spm_Tcdf(t,df(2));
        end
    end

elseif  STAT == 'X'

    % Chi-squared: SPM(X)
    %----------------------------------------------------------------------
    for d = 1:(D + 1)
        if d > 1
            EC(d,:) = ECdF(d - 1,[df(2),Inf],t);
        else
            EC(d,:) = 1 - spm_Xcdf(t,df(2));
        end
    end

elseif  STAT == 'F'

    % SPM(F)
    %----------------------------------------------------------------------
    for d = 1:(D + 1)
        if d > 1
            EC(d,:) = ECdF(d - 1,df,t);
        else
            EC(d,:) = 1 - spm_Fcdf(t,df);
        end
    end

end

function e = ECdF(d,df,t)
% Euler characteristic density for SPM(F)
%--------------------------------------------------------------------------
p = min(df(1),256);
m = min(df(2),256);
e = (log(2)/pi)^(d/2)*2*factorial(d-1)*gamma((p+m-d)/2)/...
    (m^((p-d)/2)*gamma(p/2)*gamma(m/2))*(1+(p*t)/m).^(-(p+m-2)/2);

% polynomial
%--------------------------------------------------------------------------
Si    = 0;
for i = 0:(d-1)
    Sj    = 0;
    for j = 0:min(i,d-1-i)
        b = (p+m-d)/2 + j - 1;
       Sj = Sj + K(b,j)*K(m-1,i-j)*K(p-1,d-1-i-j)*(m^(-i));
    end
    Si = Si + Sj*((-1).^(d-1-i))*((p*t).^(i+(p-d)/2));

end
e = e.*Si;

return

function g = K(b,a)
% polynomial coefficients
%--------------------------------------------------------------------------
if b - a + 1 > 0
    g = gammaln(b+1) - (gammaln(a+1) + gammaln(b-a+1));
    g = exp(g);
else
    g = 0;
end

return

% NB: numerical test
%==========================================================================
t    = linspace(1,8,32);
df   = [3,64];
STAT = 'T';
D    = 3;
EC   =  spm_ECdensity(STAT,t,df);
ec   = spm_EC_density(STAT,t,df,3);

subplot(2,2,1), plot(t,EC',t,ec',':k'),axis square, title('Old and new')
subplot(2,2,2), plot(t,ec'),axis square, title('New')

