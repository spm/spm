function f0 = push(f1, phi, dm0, sett)
% GPU single precision push
% FORMAT f0 = push(f1, phi, dm0, sett)
% f1   - 3D float array
% phi  - 4D float array (dim(4)=3)
% dm0  - Output dimensions
% sett - Settings
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


persistent kernel

if nargin<4, sett = pp_settings; end

if true
    usize = @uint64;
else
    usize = @uint32;
end

% Input image dimensions
dm1  = [size(f1) 1 1];
if prod(dm1(4:end))>1, error('Can only handle a single image.'); end
dm1  = usize(dm1(1:3));
N    = usize(prod(dm1));

% 3rd argument: Output dimensions
if nargin<3
    dm0 = dm1;
else
    dm0 = [dm0 1 1 1];
    if prod(dm0(4:end))~=1, error('Can only handle a single image.'); end
    if any(dm0~=round(dm0)) || any(dm0<0), error('Inappropriate output dimensions.'); end
    dm0 = usize(dm0(1:3));
end

% More input dimensions
dm = size(phi);
if numel(dm)~=4, error('Inappropriate dimensions for phi.'); end
if dm(4)~=3,     error('Incorrect dim(4) for phi.');         end
if ~all(dm1==dm(1:3)), error('Dimension mismatch.');    end


if isa(f1,'gpuArray') || isa(phi,'gpuArray')
    if isempty(kernel)
        % Load the kernel
        ptxfile = ptxlocation('pushpull'); % ptx file
        cproto  = 'float *, const float *, const float *'; % C prototype
        funchan = 'push_element';                          % Function handle
        kernel  = parallel.gpu.CUDAKernel(ptxfile, cproto, funchan);
        kernel.ThreadBlockSize = kernel.MaxThreadsPerBlock;
    end

    kernel.GridSize = [ceil(N/prod(kernel.ThreadBlockSize)),1];
    f0 = zeros(dm0,'like',f1);
    setConstantMemory(kernel, 'n1', usize(N), 'd0',usize(dm0), ...
                      'bnd',int32(sett.bnd), 'dp',usize(sett.deg+1), 'ext',int32(sett.ext));
    f0 = feval(kernel, f0, phi, f1);
else
    loadlib('pushpull');
    f0 = zeros(dm0,'single');
    f0 = calllib('pushpull','push', f0, phi, f1, ...
                 usize(dm0), usize(N), int32(sett.bnd), usize(sett.deg(1:3)+1), int32(sett.ext));
    f0 = reshape(f0,dm0);
end
