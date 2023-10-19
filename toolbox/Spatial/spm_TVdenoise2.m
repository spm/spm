function y = spm_TVdenoise2(x, lambda, nit, y)
%
% FORMAT y = spm_TVdenoise2(x, lambda, nit, y)
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


if nargin < 3, nit = [1000 10]; end
if nargin < 2, lambda = single([1 1 1 1 1 1 1 1]*0.01); end
if nargin < 4, y = x; end

if numel(lambda)==1
    lambda = lambda*ones(1,8,'single');
end

d    = uint64([size(x) 1 1]); % Gradient dimensions
d    = d(1:3);
if numel(nit)<2
    nit = [0 nit];
end
if isa(x,'gpuArray')
    ptxfile = ptxlocation('TVdenoise2d');
    cproto  = 'float *, const float *';

    for code=1:2
        if nit(code)==0, continue; end
        % Load the kernel
        if code==1
            funchan = 'TVdenoise2d_fast';
        else
            funchan = 'TVdenoise2d';
        end
        kernel  = parallel.gpu.CUDAKernel(ptxfile, cproto, funchan);
        kernel  = threadblocks(kernel,[ceil((d(1)-2)/2), ceil((d(2)-2)/2)]);

        setConstantMemory(kernel,'d',d,'lambda',lambda.^2);
        for it=1:nit(code)
            for oj=0:2
                for oi=0:2
                    setConstantMemory(kernel,'o',uint64([oi oj]));
                    y  = feval(kernel, y, x);
                end
            end
        end
    end
else
    warning('Done nothing.');
end
