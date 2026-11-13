function B = spm_thin(B,M)
% Erosion to find medial planes
% FORMAT B = spm_thin(B,M)
% B - 3D binary array to thin (on input)
% M - Optional additional input constraining some voxels to be 1
% B - Thinned array (on output)
%
% This code is not beautiful, but it seems to work a bit. If it is deemed
% to be sufficiently useful, then it could be speeded up by recoding parts
% of it in C.
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 1994-2025 Wellcome Centre for Human Neuroimaging

D = size(B);
F = ones(3,3,3);
F(2,2,2) = 0;
r = 0:2;
w = 2.^[0:12 0 13:25]; % For indexing the lookup table
w(14) = 0;
lkp = get_lookup;
%load topol_lookup.mat

for it=1:min(size(B)) % Iterate

    deleted = 0;            % Number of "deleted" voxels
    C0 = convn(B,F,'same'); % Count the number of neighbours

    % Work over thresholds, considering furthest voxels first. i.e.
    % voxels with fewer than 10 neighbours, then those with fewer
    % than 11, fewer than 12, etc, up to those with fewer than 26
    % neighbours. If a voxel has all 26 neighbours, then it will be
    % in the middle so not a good candidate for erosion.
    for thresh = 9:25

        C = B & (C0<=thresh);                    % Candidate voxels
        C = C(2:(D(1)-1),2:(D(2)-1),2:(D(3)-1)); % Simplify by not dealing with the edges
        ind = find(C(:))';                       % Find the candidates

        for subit=1:32 % Iterate until no more voxels can be removed
            ind = fliplr(ind); % Reverse the directions at each sub-iteration
            cnt = 0;           % Count of removed voxels
            for ii=ind         % Loop over candidates
                [i,j,k] = ind2sub(D-2,ii); % Get indices of current candidate
                b = B(i+1,j+1,k+1);        % Value of current voxel
                if nargin>=2, m = M(i+1,j+1,k+1); else, m = false; end

                if b && ~m
                    c  = B(i+r,j+r,k+r); % Get 3x3x3 patch
                    sc = sum(c(:))-1;    % Count number of neighbours
                    if sc>8 || sc<5
                        % Not entirely sure of the best values to use regarding the
                        % number of neighbours to use. A complete surface should have
                        % eight neighbours, whereas I think the edge should have about
                        % five. If a point has fewer than five neighbours, then I think
                        % we can assume that it is not part or a surface.
                        %
                        % Mid-surface   Edge
                        %    * * *      . . .
                        %    * o *      * o *
                        %    * * *      * * *

                        % Convert the pattern of 0s and 1s to an integer and read the
                        % lookup table at this value
                        d = lkp(w*c(:)+1);
                        if ~d
                            % If the lookup table says flipping from 1 to 0 leaves the
                            % topology unchanged, then set to 0
                            cnt = cnt+1;
                            B(i+1,j+1,k+1) = 0;
                            deleted = deleted + 1;
                        end
                    end
                end
            end
            if cnt==0, break; end
        end
        fprintf('%3d', subit-1);
    end
    fprintf('\n');
    if deleted==0, break; end
end
fprintf('\n');
end


function lkp = get_lookup
    % Read the lookup tabe, which is stored as a series of bytes
    % to save disk space. This is converted to a vector of 2^26
    % binary values (2^26 allows all configurations of 26 neighbours).
    [dr,~,~] = fileparts(mfilename);
    load(fullfile(dr,'topology.mat'),'numbers');
    lkp      = zeros(2^26,1,'logical');
    for bit=1:8
        lkp(bit:8:end) = bitget(numbers,bit,'uint8');
    end
end

