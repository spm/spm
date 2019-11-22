function matreg = Tikhonov_rank_def(mat, rank, lambda)

% Apply Tikhonov regularisation to rank-deficient matrix

% mat: square matrix

% rank: number of singular values to be considered in inversion

% lambda: regularisation parameters

% OH, Sep 2018


% SVD of input matrix

[U,S,V] = svd(mat);


% get singular values

s = diag(S);


% take only relevant values

s2 = s(1:rank);


lambda = (lambda/100) * (sum(s2)/rank);

% regularise eigenvalues with Tikhonov and invert

s2 = s2 ./ (s2.^2 + lambda);


% reconstitute regularised inverse matrix

matreg = V(:,1:rank)*diag(s2)*U(:,1:rank)';