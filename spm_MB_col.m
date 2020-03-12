function [col,bol,msz] = spm_MB_col(n)
% FORMAT [col,bol,msz] = spm_MB_col(n)
% Return colours and marker size for number of partitions
% n  - number of partitions
%__________________________________________________________________________
% Copyright (C) 2019-2020 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MB_col.m 7799 2020-03-12 17:23:14Z karl $

% Marker colour and size
%--------------------------------------------------------------------------
s = rand('twister');
rand('twister',1);
msz   = fix(16 + 64/n);
for k = 1:max(n,1)
    bol{k} = spm_softmax(log(rand(3,1))*2);
    col{k} = bol{k}*(1 - 1/2) + ones(3,1)/2;
end
rand('twister',s);
