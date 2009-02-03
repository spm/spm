function r = spm_gamrnd(a,b,varargin)
% Random arrays from gamma distribution
% FUNCTION r = spm_gamrnd(a,b,m,n,...)
%
% a       - shape parameter
% b       - scale parameter
% m,n,... - dimensions of the output array [optional]
%
% r       - array of random numbers chosen from the gamma distribution
%__________________________________________________________________________
%
% Reference
% 
% George Marsaglia and Wai Wan Tsang, "A Simple Method for Generating Gamma
% Variables": ACM Transactions on Mathematical Software, Vol. 26, No. 3,
% September 2000, Pages 363-372
%
% Part of RANDRAW - Efficient Random Variates Generator:
% http://www.mathworks.com/matlabcentral/fileexchange/7309
% Alex Bar Guy  &  Alexander Podgaetsky
%
% Version 1.1
%
% Any comments and suggestions please send to:
%    alex@wavion.co.il
% 
% These programs are distributed in the hope that they will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Alex Bar Guy & Alexander Podgaetsky
% $Id: spm_gamrnd.m 2685 2009-02-03 19:16:00Z guillaume $

error(nargchk(2,Inf,nargin));

if nargin == 2
    sampleSize = [1 1];
else
    sampleSize = [varargin{1:end}];
end

if a < 1
    % If a<1, one can use GAMMA(a)=GAMMA(1+a)*UNIFORM(0,1)^(1/a);
    r = spm_gamrnd(1+a, 1, sampleSize) .* (rand(sampleSize).^(1/a));

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

end % if a < 1, else ...

r = b .* r;
