function z = img_t2z(t,df, tails)
% this function was adapted from:
% https://github.com/jmcarp/fmri-pipe/blob/master/t2z.m

if nargin < 3
    tails = 2;
end

z = zeros(size(t));

% dfx = ones(size(t)) * df;
logp = t2logp(t, df);

isbig = logp < -14.05 & imag(logp) == 0;
issmall = ~isbig;

logp(issmall) = log(1 - tcdf(t(issmall), df));

if tails == 2
    logp = logp + log(2);
end

z(isbig) = logp2zbig(logp(isbig));
z(issmall) = logp2zsmall(logp(issmall));

function logp = t2logp(t, df)

lbeta = betaln(1 / 2, df ./ 2);
logp = log( ((3 .* df .^ 2 ./ ( (df + 2) .* (df + 4) .* t .^ 2) - df ./ (df + 2)) ./ t .^ 2 + 1) ./ (sqrt(df) .* t)) ...
    - ((df - 1) ./ 2) .* log(1 + t .^ 2 ./ df) - lbeta;

function z = logp2zsmall(logp)

p = exp(logp);

p(p > 0.5) = 0.5;
z = norminv(1 - p);
z(isnan(z)) = 0;

function z = logp2zbig(logp)

z = sqrt(-2 * logp - log(2 * pi));

for zidx = 1 : 3

    z = sqrt(-2 * logp - log(2 * pi) - 2 * log(z) + ...
        2 * log(1 - z .^ -2 + 3 * z .^ -4));

end

