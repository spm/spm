function y = spm_TVdenoise(x, vox, lambdap, lambdal, nit, y)
% Joint total variation denoising of 3D volumes
% FORMAT y = spm_TVdenoise(x, vox, lambdap, lambdal, nit, y)
% x        - a 3D or 4D array/gpuArray of floating point data
% vox      - voxel sizes [1 1 1]
% lambdap  - regularisation of each channel (along 4th dimension) [1]
% lambdal  - reciprocals of variances (along 4th dimension) [1]
% nit      - number of iterations [100]
% y        - starting estimates [x]
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


if nargin < 2, vox     = single([1,1,1]); end
if nargin < 3, lambdap = ones(1,20,'single'); end
if nargin < 4, lambdal = ones(1,20,'single'); end
if nargin < 5, nit     = 100; end
if nargin < 6, y       = x; end

if numel(lambdap)==1, lambdap = lambdap*ones(1,8,'single'); end
if numel(lambdal)==1, lambdal = lambdal*ones(1,8,'single'); end

if true
    usize = @uint64;
else
    usize = @uint32;
end

d = usize([size(x) 1 1]); % Gradient dimensions
d = d(1:4);
if d(4)>8, error('Too many volumes.'); end

if isa(x,'gpuArray')

    % Load the kernel
    ptxfile = ptxlocation('TVdenoise3d');
    cproto  = 'float *, const float *';
    %funchan = 'TVdenoise3d_fast';
    funchan = 'TVdenoise3d';
    kernel  = parallel.gpu.CUDAKernel(ptxfile, cproto, funchan);
    kernel  = threadblocks(kernel,ceil((d(1:3)-2)/2));

    setConstantMemory(kernel,'d',usize(d), 'vox',single(vox), 'lambdap',single(lambdap.^2), 'lambdal',single(lambdal));
    for it=1:nit
        for ok=0:2
            for oj=0:2
                for oi=0:2
                    setConstantMemory(kernel,'o',usize([oi oj ok]));
                    y  = feval(kernel, y, x);
                end
            end
        end
        %fprintf('.');
    end
    %fprintf('\n');
else
    loadlib('TVdenoise3d');
    y = calllib('TVdenoise3d','TVdenoise3d', y, x, usize(d), single(vox), single(lambdap.^2), single(lambdal), uint32(nit));
    y = reshape(y,d); % Needed because callib seems to return 2D arrays only
end
