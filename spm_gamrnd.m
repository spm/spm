function r = spm_gamrnd(a,b,varargin)
% THE GAMMA DISTRIBUTION
%
% The standard form of the GAMMA distribution:
%
%   pdf(y) = y^(a-1)*exp(-y)/gamma(a);  y>=0, a>0
%   cdf(y) = gammainc(y, a);
%
%   Mean = a;
%   Variance = a;
%   Skewness = 2/sqrt(a);
%   Kurtosis = 6/a;
%   Mode = a-1;
%
% The general form of the GAMMA distribution:
%
%   pdf(y) = ((y-m)/b).^(a-1) .* exp(-(y-m)/b)/ (b*gamma(a));  y>=m; a>0; b>0
%   cdf(y) = gammainc((y-m)/b, a);  y>=m; a>0; b>0
%
%   Mean = m + a*b;
%   Variance = a*b^2;
%   Skewness = 2/sqrt(a);
%   Kurtosis = 6/a;
%   Mode = m + b*(a-1);
%
% PARAMETERS:
%   m - location
%   b - scale; b>0
%   a - shape; a>0
%
% SUPPORT:
%   y,   y>=0   - standard GAMMA distribution
%    or
%   y,   y>=m   - generalized GAMMA distribution
%
% CLASS:
%   Continuous skewed distributions
%
% NOTES:
% 1. The GAMMA distribution approaches a NORMAL distribution as a goes to Inf
% 5. GAMMA(m, b, a), where a is an integer, is the Erlang distribution.
% 6. GAMMA(m, b, 1) is the Exponential distribution.
% 7. GAMMA(0, 2, nu/2) is the Chi-square distribution with nu degrees of freedom.
%
% USAGE:
%   randraw('gamma', a, sampleSize) - generate sampleSize number
%         of variates from standard GAMMA distribution with shape parameter 'a';
%   randraw('gamma', [m, b, a], sampleSize) - generate sampleSize number
%         of variates from generalized GAMMA distribution
%         with location parameter 'm', scale parameter 'b' and shape parameter 'a';
%   randraw('gamma') - help for GAMMA distribution;
%
% EXAMPLES:
%  1.   y = randraw('gamma', [2], [1 1e5]);
%  2.   y = randraw('gamma', [0 10 2], 1, 1e5);
%  3.   y = randraw('gamma', [3], 1e5 );
%  4.   y = randraw('gamma', [1/3], 1e5 );
%  5.   y = randraw('gamma', [1 3 2], [1e5 1] );
%  6.   randraw('gamma');
%
% END gamma HELP END gama HELP
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
%
% $Id: spm_gamrnd.m 2635 2009-01-21 09:43:47Z maria $

% Method:
%
% Reference:
% George Marsaglia and Wai Wan Tsang, "A Simple Method for Generating Gamma
%   Variables": ACM Transactions on Mathematical Software, Vol. 26, No. 3,
%   September 2000, Pages 363-372

%   EFFICIENT RANDOM VARIATES GENERATOR
%
%   http://www.mathworks.com/matlabcentral/fileexchange/7309
%
%  Version 1.1 - April 2005 -  Bug fix:   Generation from binomial distribution using only 'binomial'
%                                   usage string was changed to 'binom' ( 'binomial' works too ).
%  Version 1.0 - March 2005 -  Initial version
%  Alex Bar Guy  &  Alexander Podgaetsky
%    alex@wavion.co.il

% These programs are distributed in the hope that they will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

% Any comments and suggestions please send to:
%    alex@wavion.co.il

error(nargchk(2,Inf,nargin));

if nargin == 2
    sampleSize = [1 1];
else
    sampleSize = [varargin{1:end}];
end

m = 0;

if a < 1
    % If a<1, one can use GAMMA(a)=GAMMA(1+a)*UNIFORM(0,1)^(1/a);
    r = m + b*(spm_gamrnd(1+a, 1, sampleSize)).*(rand(sampleSize).^(1/a));

else

    d = a - 1/3;
    c = 1/sqrt(9*d);

    x = randn( sampleSize );
    v = 1+c*x;

    indxs = find(v <= 0);
    while ~isempty(indxs)
        indxsSize = size( indxs );
        xNew = randn( indxsSize );
        vNew = a+c*xNew;

        l = (vNew > 0);
        v( indxs( l ) ) = vNew(l);
        x( indxs( l ) ) = xNew(l);
        indxs = indxs( ~l );
    end

    u = rand( sampleSize );
    v = v.^3;
    x2 = x.^2;
    r = d*v;

    indxs = find( (u>=1-0.0331*x2.^2) & (log(u)>=0.5*x2+d*(1-v+log(v))) );
    while ~isempty(indxs)
        indxsSize = size( indxs );

        x = randn( indxsSize );
        v = 1+c*x;
        indxs1 = find(v <= 0);
        while ~isempty(indxs1)
            indxsSize1 = size( indxs1 );
            xNew = randn( indxsSize1 );
            vNew = a+c*xNew;

            l1 = (vNew > 0);
            v( indxs1(l1) ) = vNew(l1);
            x( indxs1(l1) ) = xNew(l1);
            indxs1 = indxs1( ~l1 );
        end

        u = rand( indxsSize );
        v = v .* v .* v;
        x2 = x.*x;

        l = (u<1-0.0331*x2.*x2) | (log(u)<0.5*x2+d*(1-v+log(v)));
        r( indxs( l ) ) = d*v(l);
        indxs = indxs( ~l );
    end % while ~isempty(indxs)

    r = m + b*r;

end % if a < 1, else ...
