function s = getthreads(kernel,d)
% Size of block of threads on a CUDA kernel
% FORMAT s = getthreads(kernel,d)
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


nmax = kernel.MaxThreadsPerBlock;
d    = d(1:min(3,numel(d)));
d    = double(d);
s    = ones(size(d));
c    = cumprod(d);
n    = numel(d);
if n>1
    if c(1)>nmax
        s(1) = nmax;
    else
        s(1) = d(1);
        if n>2
            if c(2)>nmax
                s(2) = floor(nmax/c(1));
            else
                s(2) = d(2);
                if n>3
                    if c(3)>nmax
                        s(3) = floor(nmax/(c(2)));
                    else
                        s(3) = d(3);
                    end
                end 
            end
        end
    end
end
