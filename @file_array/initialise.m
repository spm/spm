function initialise(fa,v)
% Initialise file on disk
%
% This creates a file on disk with the appropriate size by explicitly
% writing data to prevent a sparse file.
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

%
% $Id: initialise.m 5432 2013-04-22 13:02:45Z guillaume $

if nargin==1, v = 0; end

%fa = subsasgn(fa, substruct('()',num2cell(size(fa))), v);

bs = 2^20;
n  = prod(size(fa)); %#ok<PSIZE>
fa = reshape(fa,n,1);
for i=1:ceil(n/bs)
    ii = ((((i-1)*bs)+1):min((i*bs),n))';
    fa = subsasgn(fa,struct('type','()','subs',{{ii}}),repmat(v,numel(ii),1));
end
