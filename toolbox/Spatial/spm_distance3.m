function d = spm_distance3(p,vx,maxit,d)
% Compute Euclidean distances from a tissue map
% FORMAT d = spm_distance3(p,vx,maxit,d)
% p     - Tissue probability map
% vx    - Voxel sizes ([1 1 1])
% maxit - Maximum number of iterations (20)
% d     - Distances.
%         - Initial estimates on input (zeros(size(p)))
%         - Output estimate on output
%
% This implementation is very slow. Vectorising it would be memory
% hungry, so (if it turns out to be useful) it should be written in C.
%
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 1994-2025 Wellcome Centre for Human Neuroimaging

if nargin < 4, d     = zeros(size(p)); end
if nargin < 3, maxit = 20;             end
if nargin < 2, vx    = [1 1 1];        end

% Determine distances (dp) of neighbouring voxels,
% as well as their indices (ind).
[dx,dy,dz] = ndgrid((-1:1)*vx(1),(-1:1)*vx(2),(-1:1)*vx(3));
dp  = sqrt(dx.^2+dy.^2+dz.^2);
ind = [1:13 15:26]';
dp  = dp(ind);

for it=1:maxit % Iterate
    changes = false;

    for i1=1:3
        for j1=1:3
            for k1=1:3

                for i=(1+i1):3:(size(p,1)-1)
                    for j=(1+j1):3:(size(p,2)-1)
                        for k=(1+k1):3:(size(p,3)-1)
                            % 3x3x3 patch of distances around voxel i,j,k
                            patch    = d((-1:1)+i,(-1:1)+j,(-1:1)+k);

                            % Compute distance from tissues
                            s        = (1-p(i,j,k)).*min(patch(ind) + dp);

                            % Check for any changes
                            if s~= d(i,j,k), changes = true; end

                            % Update the distance
                            d(i,j,k) = s;
                        end
                    end
                end
            end
        end
    end
    fprintf('.');
    if ~changes, break; end
end
fprintf('\n');

