function f = lbessi(nu,z)
% GPU single precision f = log(besseli(nu, z))
persistent kernel
if isempty(kernel)
    % Compilation needs the cuda compiler (nvcc) to be installed.
    % On architectures I tried:
    %
    % Linux (after adding /usr/local/cuda-10.0/bin/ to the path):
    %   nvcc -ptx lbessi.cu -o lbessi.ptxa64
    %
    % Windows (cygwin with /cygdrive/c/Program Files (x86)/ ...
    %          Microsoft Visual Studio/2019/Community/VC/ ...
    %          Tools/MSVC/14.29.30133/bin/Hostx86/x64
    %          and  /cygdrive/c/Program Files/ ...
    %          NVIDIA GPU Computing Toolkit/CUDA/v10.1/bin/nvcc
    %          added to the path):
    %  nvcc -ptx lbessi.cu -o lbessi.ptxw64
    %

    % Load the kernel
    ptxfile = ptxlocation('lbessi');
    disp(ptxfile)
    dir(ptxfile)
    cproto  = 'float *, const float, const float *, const unsigned int';
    funchan = 'lbessi_element';
    kernel = parallel.gpu.CUDAKernel(ptxfile, cproto, funchan);
    kernel.ThreadBlockSize = [kernel.MaxThreadsPerBlock,1,1];
end
N = numel(z);
kernel.GridSize = [ceil(N/kernel.MaxThreadsPerBlock),1];

f = zeros(size(z),'like',z);
f = feval(kernel, f, nu, z, N);

