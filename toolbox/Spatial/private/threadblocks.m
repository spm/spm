function kernel = threadblocks(kernel,d)
% Set the size of a block of threads and grid on a CUDA kernel
% FORMAT kernel = threadblocks(kernel,d)
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging

threads = getthreads(kernel,d);
blocks  = ceil(double(d)./threads);
kernel.ThreadBlockSize = threads;
kernel.GridSize = blocks;
