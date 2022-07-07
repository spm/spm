function [col,bol,msz] = spm_MB_col(n)
% Return colours and marker size for number of partitions
% FORMAT [col,bol,msz] = spm_MB_col(n)
% n  - number of partitions
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2019-2022 Wellcome Centre for Human Neuroimaging


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
