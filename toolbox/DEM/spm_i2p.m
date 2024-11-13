function P = spm_i2p(X,s)
% Converts an array of indices into a cell array of probabilities
% FORMAT P = spm_i2p(X,s)
% X   - (m,n) : 0 < X(i,j) < s
% P   - {m,n} : P{i,j}(X(i,j),j) = 1
% 
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2022-2023 Wellcome Centre for Human Neuroimaging

% create logical cell array of one-hot vectors
%--------------------------------------------------------------------------
[m,n]  = size(X);
P      = cell(m,n);
[P{:}] = deal(false(s,1));
for i = 1:m
    for j = 1:n
        P{i,j}(X(i,j)) = true;
    end
end

return


