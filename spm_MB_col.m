function [col,bol,msz] = spm_MB_col(n)
% FORMAT [col,bol,msz] = spm_MB_col(n)
% Return colours and marker size for number of partitions
% n  - number of partitions
%__________________________________________________________________________
% Copyright (C) 2019-2020 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MB_col.m 7768 2020-01-07 11:37:19Z spm $

% Marker colour and size
%--------------------------------------------------------------------------
s = rand('twister');
rand('twister',1);
msz   = fix(16 + 64/n);
for k = 1:n
    bol{k} = spm_softmax(log(rand(3,1))*2);
    col{k} = bol{k}*(1 - 1/2) + ones(3,1)/2;
end
rand('twister',s);
