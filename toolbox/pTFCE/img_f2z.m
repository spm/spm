function z = img_f2z(f, d1, d2, tails)

if d1 == 1
    z = img_t2z(sqrt(f), d2, tails);
    return
end

if nargin < 4
    tails = 2;
end

z = zeros(size(f));

logp = f2logp(f, d1, d2);

isbig = logp < -14.05 & imag(logp) == 0;
issmall = ~isbig;

logp(issmall) = log(1 - fcdf(f(issmall), d1, d2));

if tails == 1
    logp = logp - log(2);
end

z(isbig) = logp2zbig(logp(isbig));
z(issmall) = logp2zsmall(logp(issmall));

function logp = f2logp(f, d1, d2)

a = d1 / d2;
m = (d1 + d2) / 2;
n = 1 - d1 / 2;

loggam = (d1 / 2) * log(d1 / d2) - betaln(d2 / 2, d1 / 2);

top = 1;
bot = n + m - 1;
iter = 0;

for iteridx = 1 : 20
    iter = iter + top * (f .^ (-(n + iteridx - 1)) / (a ^ (iteridx) * bot));
    top = top * (n - 1 + iteridx) * -1;
    bot = bot * (n + m - 1 + iteridx);
end

logp = loggam - (m - 1) * log(1 + a * f) + log(iter);
if iter <= 0
    logp = inf;
end

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

