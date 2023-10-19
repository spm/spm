function f = lbessi(nu,z)
% GPU single precision f = log(besseli(nu, z))
% FORMAT f = lbessi(nu,z)
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


persistent kernel
if isempty(kernel)
    % Load the kernel
    ptxfile = ptxlocation('lbessi');
    cproto  = 'float *, const float, const float *, const unsigned int';
    funchan = 'lbessi_element';
    kernel  = parallel.gpu.CUDAKernel(ptxfile, cproto, funchan);
    kernel.ThreadBlockSize = [kernel.MaxThreadsPerBlock,1,1];
end
N = numel(z);
kernel.GridSize = [ceil(N/kernel.MaxThreadsPerBlock),1];

f = zeros(size(z),'like',z);
f = feval(kernel, f, nu, z, N);
