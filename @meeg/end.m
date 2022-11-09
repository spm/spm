function en = end(this,k,n)
% Overloaded end function for meeg objects.
%__________________________________________________________________________

% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


dim = size(this.data);
if k>length(dim)
    en = 1;
else
    if n<length(dim)
    dim = [dim(1:(n-1)) prod(dim(n:end))];
    end
    en = dim(k);
end
