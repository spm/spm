function p = spm_Q_perm(Q)
% Return a cell of permutation indices for separating matrices
% FORMAT p = spm_Q_perm(Q)
%__________________________________________________________________________
 
% John Ashburner & Karl Friston
% Copyright (C) 2009-2022 Wellcome Centre for Human Neuroimaging
 
 
% add cell arrays of components
%--------------------------------------------------------------------------
A   = 0;
if iscell(Q)
    for i = 1:length(Q)
        A = A + ~~Q{i};
    end
else
    A = ~~Q;
end
 
% find partition
%--------------------------------------------------------------------------
A = A + A';
p = {};
i = find(any(A,1),1);
while ~isempty(i)
    k     = [];
    while length(k) < length(i)
        k = i;
        j = find(any(A(:,i),2));
        i = j;
    end
    p{end + 1} = i;
    A(i,i)     = 0;
    i          = find(any(A,1),1);
end
