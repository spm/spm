function lkp = spm_topo_lookup
% Create a topology change lookup table
% FORMAT lkp = spm_topo_lookup
% lkp - The lookup table of 1s and 0s.
%
% The table is used for determining whether erosion can be
% done without changing local topology. In addition to returning
% lkp, it is also saved in a topol_lookup.mat file.
%
% The code is not very efficient, so the lookup table takes about
% a day and a half to compute.
%__________________________________________________________________________

% John Ashburner
% Copyright (C) 1994-2025 Wellcome Centre for Human Neuroimaging

lkp = lookup;
save topo_lookup.mat lkp
end

function lkp = lookup
    % Loop through all 2^26 possible patterns of neighbours and
    % return a vector of ones and zeros that indicate whether
    % flipping the voxel at the centre of a 3x3x3 patch will change
    % the Euler characteristic of that patch. The EC indicates how
    % many objects there would be, so when doing skeletonisation, we
    % do not want to change this. There are 2^26 possible arrangements
    % of 26 neighbours, so the lookup table has that length.
    N   = 2^26;
    lkp = zeros([N,1],   'logical');
    X   = zeros([3,3,3], 'logical');
    ind = [1:13 15:27];
    for i=1:N
        X(ind)  = bitget(i-1,1:26,'Int32');
        lkp(i)  = topo_change(X);

        % Show a few dots so we know it is working.
        if ~rem(i,262144*2),  fprintf('.'); end
        if ~rem(i,8388608*4), fprintf('\n'); end
    end
    fprintf('\n');
end


function t = topo_change(X)
    % Does flipping the centre voxel of a 3x3x3 cube change the
    % topology?

    % Compute Euler characteristic (divided by 2) if
    % the centre voxel has a value of 0.
    X(14)   = false;
    [v,e,f,p] = VEFP_KW(X);
    ec0     = v - e + f - p;

    % Compute Euler characteristic if
    % the centre voxel has a value of 1.
    X(14)   = true;
    [v,e,f,p] = VEFP_KW(X);
    ec1     = v - e + f - p;

    % Return true if the two ECs differ (false if flipping
    % the voxel does not change the EC).
    t       = (ec0~=ec1);
end

function [v,e,f,p] = VEFP_KW(X0)
    % Compute statistics to derive the Euler characteristic
    % of a 3x3x3 cube of voxels, which is done as described in
    %     Worsley, K.J., 1996. The geometry of random images.
    %     Chance, 9(1), pp.27-40.
    %     https://static.wikitide.net/biuwiki/9/9a/1996_Worsley.pdf
    %
    % See also src/spm_resels_vol.c
    %
    % The outputs are:
    % v - number of vertices
    % e - number of edges
    % f - number of faces
    % p - number of polyhedra

    X = zeros(5,5,5); % Pad to simplify things
    X(2:4,2:4,2:4) = X0;
    p = 0;
    f = 0;
    e = 0;
    v = 0;
    for i=1:4
        ir = i+(0:1);
        for j=1:4
            jr = j+(0:1);
            for k=1:4
                kr = k+(0:1);

                % Polyhedra
                p  = p + X(i,j,k);

                % Faces
                f  = f + inside(X(ir,j ,k ));
                f  = f + inside(X(i ,jr,k ));
                f  = f + inside(X(i ,j, kr));

                % Edges
                e  = e + inside(X(i ,jr,kr));
                e  = e + inside(X(ir,j ,kr));
                e  = e + inside(X(ir,jr,k ));

                % Vertices
                v  = v + inside(X(ir,jr,kr));
            end
        end
    end
end

function t = inside(x)
    t = any(x(:));
end

