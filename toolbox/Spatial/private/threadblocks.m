function kernel = threadblocks(kernel,d)
threads = getthreads(kernel,d);
blocks  = ceil(double(d)./threads);
kernel.ThreadBlockSize = threads;
kernel.GridSize = blocks;

