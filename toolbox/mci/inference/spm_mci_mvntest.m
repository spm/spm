function [stats] = spm_mci_mvntest (X,df)
% Test for multivariate normality
% FORMAT [stats] = spm_mci_mvntest (X,df)
%
% X         [N x d] data matrix containing N IID samples of d-variate data
% df        Degrees of freedom e.g. df <= N
%
% stats
% .p        p-value for multivariate normality
%           e.g. with p < 0.05 one can reject null hypothesis of normality 
% .W(j)     Shapiro-Wilks statistic for jth variate
% .Z(j)     Equivalent standardised Gaussian variate for W(j)
% .pusw(j)  p-value for normality of jth variate (Shapiro-Wilks)
% .puks(k)  p-value for normality of jth variate (Kolmogorov-Smirnoff)
% .k        kurtosis of jth variate
% .s        skewness of jth variate
%
% This function is adapted from the matlab function Roystest.m:
%
% Trujillo-Ortiz, A., R. Hernandez-Walls, K. Barba-Rojo and
%   L. Cupul-Magana. (2007). Roystest:Royston's Multivariate Normality Test.   
%   A MATLAB file. [WWW document]. URL http://www.mathworks.com/
%   matlabcentral/fileexchange/loadFile.do?objectId=17811
%
% Royston, J.P. (1992). Approximating the Shapiro-Wilk W-Test for non-
%      normality. Statistics and Computing, 2:117-119.
%      121-133.
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_mci_mvntest.m 7679 2019-10-24 15:54:07Z spm $

alpha = 0.05;

if nargin < 1,
    error('Requires at least one input argument.');
    return;
end

if nargin < 2 | isempty(df)
    [n,p] = size(X);
else
    n = df;
    p = size(X,2);
end

if (n <= 3),
    error('n is too small.');
    return,
elseif (n >= 4) && (n <=11),
    x = n;
    g = -2.273 + 0.459*x;
    m = 0.5440 - 0.39978*x + 0.025054*x^2 - 0.0006714*x^3;
    s = exp(1.3822 - 0.77857*x + 0.062767*x^2 - 0.0020322*x^3); 
    for j = 1:p,
        W(j) = ShaWilstat(X(:,j));
        Z(j) = (-log(g - (log(1 - W(j)))) - m)/s;
    end
elseif (n >= 12) && (n <=2000),
    x = log(n);
    g = 0;
    m = -1.5861 - 0.31082*x - 0.083751*x^2 + 0.0038915*x^3;
    s = exp(-0.4803 -0.082676*x + 0.0030302*x^2);  
    for j = 1:p,
        W(j) = ShaWilstat(X(:,j));
        Z(j) = ((log(1 - W(j))) + g - m)/s;
    end
else
    error('n is not in the proper size range.'); %error('n is too large.');return,
    return,
end

for j = 1:p,
    R(j) = (norminv((normcdf( - Z(j)))/2))^2;
end

u = 0.715;
v = 0.21364 + 0.015124*(log(n))^2 - 0.0018034*(log(n))^3;
l = 5;
C = corrcoef(X); %correlation matrix
NC = (C.^l).*(1 - (u*(1 - C).^u)/v); %transformed correlation matrix
T = sum(sum(NC)) - p; %total
mC = T/(p^2 - p); %average correlation
e = p/(1 + (p - 1)*mC); %equivalent degrees of freedom
H = (e*(sum(R)))/p; %Royston's statistic
P = 1 - chi2cdf(H,e); %P-value

stats.p=P;
stats.W=W;
stats.Z=Z;
for j=1:p,
    stats.k(j)=kurtosis(X(:,j));
    stats.s(j)=skewness(X(:,j));
    stats.pusw(j)=1-spm_Ncdf(stats.Z(j),0,1);
    zj=(X(:,j)-mean(X(:,j)))/std(X(:,j));
    [h,pks]=kstest(zj);
    stats.puks(j)=pks;
end

%-----------------------------------
function [W] = ShaWilstat(x)
%SHAWILTEST Shapiro-Wilk' W statistic for assessing a sample normality.
% This m-file is generating from the Fortran Algorithm AS R94 (Royston,
% 1995) [URL address http://lib.stat.cmu.edu/apstat/181]. Here we take only
% the procedure to generate the Shapiro-Wilk's W statistic, needed to feed
% the Royston's test for multivariate normality. Here, we present both the
% options for the sample kurtosis type: 1) Shapiro-Francia for leptokurtic,
% and 2) Shapiro-Wilk for the platykurtic ones. Do not assume that the
% result of the Shapiro-Wilk test is clear evidence of normality or non-
% normality, it is just one piece of evidence that can be helpful.
%
% Input:
%      x - data vector (3 < n < 5000)
%
% Output:
%      W - Shapiro-Wilk's W statistic
%
% Reference:
% Scholz, F.W. and Stephens, M.A. (1987), K-Sample Anderson-Darling Tests.
%     Journal of the American Statistical Association, 82:918-924.
%

x  =  x(:); %put data in a column vector
n = length(x); %sample size

if n < 3,
ois   error('Sample vector must have at least 3 valid observations.');
end

if n > 5000,
    warning('Shapiro-Wilk statistic might be inaccurate due to large sample size ( > 5000).');
end

x = sort(x); %sorting of data vector in ascending order
m = norminv(((1:n)' - 3/8) / (n + 0.25));
w = zeros(n,1); %preallocating weights

if kurtosis(x) > 3, %Shapiro-Francia test is better for leptokurtic samples
    w = 1/sqrt(m'*m) * m;
    W = (w' * x) ^2 / ((x - mean(x))' * (x - mean(x)));
else %Shapiro-Wilk test is better for platykurtic samples
    c = 1/sqrt(m' * m) * m;
    u = 1/sqrt(n);
    p1 = [-2.706056,4.434685,-2.071190,-0.147981,0.221157,c(n)];
    p2 = [-3.582633,5.682633,-1.752461,-0.293762,0.042981,c(n-1)];

    w(n) = polyval(p1,u);
    w(1) = -w(n);

    if n == 3,
        w(1) = 0.707106781;
        w(n) = -w(1);
    end

    if n >= 6,
        w(n-1) = polyval(p2,u);
        w(2) = -w(n-1);
    
        ct =  3;
        phi = (m'*m - 2 * m(n)^2 - 2 * m(n-1)^2) / ...
                (1 - 2 * w(n)^2 - 2 * w(n-1)^2);
    else
        ct = 2;
        phi = (m'*m - 2 * m(n)^2) / (1 - 2 * w(n)^2);
    end

    w(ct:n-ct+1) = m(ct:n-ct+1) / sqrt(phi);

    W = (w' * x) ^2 / ((x - mean(x))' * (x - mean(x)));
end

return,
