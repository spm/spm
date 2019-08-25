function [col,bol,msz] = spm_MB_col(n)
% FORMAT [col,bol,msz] = spm_MB_col(n)
% returns colours and market size for number of partitions
% n  - number of partitions
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MB_col.m 7655 2019-08-25 20:10:20Z karl $

% Marker color and size
%--------------------------------------------------------------------------
rng(1);
msz   = fix(16 + 64/n);
for k = 1:n
    bol{k} = spm_softmax(log(rand(3,1))*2);
    col{k} = bol{k}*(1 - 1/2) + ones(3,1)/2;
end
