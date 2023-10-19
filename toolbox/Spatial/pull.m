function f1 = pull(f0, phi, sett)
% GPU single precision pull
% FORMAT f1 = pull(f0, phi, sett)
% f0   - 3D float array
% phi  - 4D float array (dim(4)=3)
% sett - Settings
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


persistent kernel

if nargin<3, sett = pp_settings; end

if true
    usize = @uint64;
else
    usize = @uint32;
end

% Input image dimensions
d0  = [size(f0) 1 1];
if prod(d0(4:end))>1, error('Can only handle a single image.'); end
d0  = usize(d0(1:3));


% Output dimensions
d1 = size(phi);
if numel(d1)~=4, error('Inappropriate dimensions for phi.'); end
if d1(4)~=3,     error('Incorrect dim(4) for phi.');         end
d1 = usize(d1(1:3));
n1 = usize(prod(d1));

if isa(f0,'gpuArray') || isa(phi,'gpuArray')
    if isempty(kernel)
        % Load the kernel
        ptxfile = ptxlocation('pushpull'); % ptx file
        cproto  = 'float *, const float *, const float *'; % C prototype
        funchan = 'pull_element';                          % Function handle
        kernel  = parallel.gpu.CUDAKernel(ptxfile, cproto, funchan);
        kernel.ThreadBlockSize = [kernel.MaxThreadsPerBlock,1,1];
       %kernel.ThreadBlockSize = kernel.MaxThreadsPerBlock;
    end

    f1 = zeros(d1,'single','gpuArray');
    kernel.GridSize = ceil(n1/kernel.MaxThreadsPerBlock);
    setConstantMemory(kernel,'d0',usize(d0), 'n1',usize(n1), ...
                      'bnd',int32(sett.bnd), 'dp',usize(sett.deg(1:3)+1), 'ext',int32(sett.ext));
    f1 = feval(kernel, f1, phi, f0);
else
    loadlib('pushpull');
    f1 = zeros(d1,'single');
    f1 = calllib('pushpull','pull', f1, phi, f0, ...
                 usize(d0), usize(n1), int32(sett.bnd), usize(sett.deg(1:3)+1), int32(sett.ext));
end
