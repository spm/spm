function [M1] = spm_eeg_inv_icp(data1,data2,fid1,fid2,Fmri,Fhsp,aff)
% Iterative Closest Point (ICP) registration algorithm.
% Surface matching computation: registration from one 3D surface onto
% another 3D surface 
%
% FORMAT [M1] = spm_eeg_inv_icp(data1,data2,[fid1],[fid2],[Fmri],[Fhsp],[aff])
% Input:
% data1      - locations of the first set of points corresponding to the
%              3D surface to register onto [3 x n]
% data2      - locations of the second set of points corresponding to the
%              second 3D surface to be registered [3 x p]
% fid1       - sMRI fiducials [default: []]
% fid2       - sens fiducials [default: []]
% Fmri       - graphics handle for sMRI points [default: none]
% Fhsp       - graphics handle for headshape [default: none]
% aff        - flag for 12-parameter affine transform [default: 0]
%
% Output:
% M1         - 4 x 4 affine transformation matrix for sensor space
%
% Landmarks (fiducials) based registration
% Fiducial coordinates must be given in the same order in both files
%
%--------------------------------------------------------------------------
% Adapted from code available at http://www.csse.uwa.edu.au/~ajmal/code/icp.m
% written by Ajmal Saeed Mian {ajmal@csse.uwa.edu.au}, Computer Science, 
% The University of Western Australia. The code may be used, modified and 
% distributed for research purposes with acknowledgement of the author and 
% inclusion this copyright information.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout & Guillaume Flandin
% $Id: spm_eeg_inv_icp.m 2345 2008-10-16 11:31:35Z guillaume $

% use figure and fiducials if specified
%--------------------------------------------------------------------------
try, fid1; catch, fid1 = []; end
try, fid2; catch, fid2 = []; end
try, aff;  catch, aff  = 0;  end


% initialise rotation and translation of sensor space
%--------------------------------------------------------------------------
M1    = speye(4,4);
tri   = delaunayn(data1');
for k = 1:64

    % find nearest neighbours
    %----------------------------------------------------------------------
    [corr, D] = dsearchn(data1', tri, data2');
    corr(:,2) = [1 : length(corr)]';
    i         = find(D > 32);
    corr(i,:) = [];
    M         = [fid1 data1(:,corr(:,1))];
    S         = [fid2 data2(:,corr(:,2))];

    % apply and accumlate affine scaling
    %----------------------------------------------------------------------
    if aff
        M     = pinv([S' ones(length(S),1)])*M';

        M     = [M'; 0 0 0 1];    
    else
        % 6-parameter affine (i.e. rigid body)
        %------------------------------------------------------------------
        M        = spm_eeg_inv_rigidreg(M,S);
    end

    data2     = M*[data2; ones(1,size(data2,2))];
    data2     = data2(1:3,:);
    fid2      = M*[fid2; ones(1,size(fid2,2))];
    fid2      = fid2(1:3,:);
    M1        = M*M1;

    % plot
    %----------------------------------------------------------------------
    try
        set(Fmri,'XData',data1(1,:),'YData',data1(2,:),'ZData',data1(3,:));
        set(Fhsp,'XData',data2(1,:),'YData',data2(2,:),'ZData',data2(3,:));
        drawnow
    end

    if (norm(M)-1)< 1e-3
        break;
    end
end


if aff == 2
    % Enforce uniform scaling for MEG case in order not to distort the head
    %----------------------------------------------------------------------
    [U, L, V] = svd(M1(1:3, 1:3));
    L = eye(3)*mean(diag(L));
    M1(1:3,1:3) =U*L*V';
end