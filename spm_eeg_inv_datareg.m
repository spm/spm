function [M1, sens, fid] = spm_eeg_inv_datareg(sens, fid, datareg, template)
% Rigid registration of the EEG/MEG data and sMRI spaces
% rigid co-registration
%           1: fiducials based (3 landmarks: nasion, left ear, right ear)
%           2: surface matching between sensor mesh and headshape
%           (starts with a type 1 registration)
%
% FORMAT [eeg2mri,sen_reg,fid_reg,hsp_reg,orient_reg,mri2eeg,hsp2eeg]
%          = spm_eeg_inv_datareg(sensors,fid_eeg,fid_mri,headshape,scalpvert,megorient,template)
% Input:
%
% sensors       - sensors struct
% fid_eeg       - the fiducial coordinates in sensor
%               space ([Nx Ny Nz ; LEx LEy LEz ; REx REy REz])
% fid_mri       - the fiducial coordinates in sMRI
%               space ([nx ny nz ; lex ley lez ; rex rey rez])
% headshape     - the headshape point coordinates in
%               sensor space ([hx1 hy1 hz1 ; hx2 hy2 hz2 ; ...])
% scalpvert     - the vertices coordinates of the scalp
%               tesselation in mri space ([vx1 vy1 vz1 ; vx2 vy2 vz2 ;
%               ...])
% template      - 0/1 switch to enable affine mri2eeg transform
%
% IMPORTANT: all the coordinates must be in the same units (mm).
%
% Output:
%
% If a template is used, the senor locations are transformed using an
% affine (rigid body) mapping.  If headshape locations are supplied
% this is generalized to a full twelve parameter affine mapping (n.b.
% this might not be appropriate for MEG data).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id: spm_eeg_inv_datareg.m 1407 2008-04-15 10:09:18Z vladimir $

% Check input arguments
%==========================================================================

% spm_eeg_inv_datareg(sensors,fid_eeg,fid_mri,headshape,scalpvert,megorient,template)
%--------------------------------------------------------------------------

if exist('template')~=1
    template  = 0;
end

fid_mri = datareg.fid_mri;
nfid  = size(fid_mri,1);

scalpvert = datareg.scalpvert;

fid.fid.pnt = fid.fid.pnt(1:nfid, :);
fid.fid.label = fid.fid.label(1:nfid);

fid_eeg = fid.fid.pnt;

% Estimate-apply rigid body transform to sensor space
%--------------------------------------------------------------------------
M1 = spm_eeg_inv_rigidreg(fid_mri', fid_eeg');

[fid_eeg, headshape, fid] = apply_transformation(M1, fid);

if isfield(sens, 'ori') % MEG

    M = headcoordinates(fid_eeg(1,:), fid_eeg(2, :), fid_eeg(3, :));
    
    M1 = inv(M) * M1;

else

    if template

        % constatined affine transform
        %--------------------------------------------------------------------------
        aff   = 1;
        for i = 1:16

            % scale
            %----------------------------------------------------------------------
            M       = pinv(fid_eeg(:))*fid_mri(:);
            M       = sparse(1:4,1:4,[M M M 1]);

            [fid_eeg, headshape, fid] = apply_transformation(M, fid);

            M1      = M*M1;

            % and move
            %----------------------------------------------------------------------
            M       = spm_eeg_inv_rigidreg(fid_mri', fid_eeg');

            [fid_eeg, headshape, fid] = apply_transformation(M, fid);

            M1      = M*M1;

        end
    else
        aff = 0;
    end
    

    % Surface matching between the scalp vertices in MRI space and
    % the headshape positions in data space
    %--------------------------------------------------------------------------
    if length(headshape)

        % load surface locations from sMRI
        %----------------------------------------------------------------------
        if size(headshape,2) > size(headshape,1)
            headshape = headshape';
        end
        if size(scalpvert,2) > size(scalpvert,1)
            scalpvert = scalpvert';
        end


        % intialise plot
        %----------------------------------------------------------------------
        h    = spm_figure('GetWin','Graphics');
        clf(h); figure(h)
        set(h,'DoubleBuffer','on','BackingStore','on');
        Fmri = plot3(scalpvert(:,1),scalpvert(:,2),scalpvert(:,3),'ro','MarkerFaceColor','r');
        hold on;
        Fhsp = plot3(headshape(:,1),headshape(:,2),headshape(:,3),'bs','MarkerFaceColor','b');
        axis off image
        drawnow

        % nearest point registration
        %----------------------------------------------------------------------
        M    = spm_eeg_inv_icp(scalpvert',headshape',fid_mri',fid_eeg',Fmri,Fhsp,aff);

        % transform headshape and eeg fiducials
        %----------------------------------------------------------------------
        [fid_eeg, headshape, fid] = apply_transformation(M, fid);
        M1        = M*M1;
    end
end

% Update the sensor locations and orientation
%--------------------------------------------------------------------------
sens = transform_sens(M1, sens);

return

%==========================================================================
function [fid_eeg, headshape, fid] = apply_transformation(M, fid)

fid = transform_headshape(M, fid);
fid_eeg = fid.fid.pnt;
headshape = fid.pnt;


%==========================================================================
function [M1] = spm_eeg_inv_icp(data1,data2,fid1,fid2,Fmri,Fhsp,aff)

% Iterative Closest Point (ICP) registration algorithm.
% Surface matching computation: registration from one 3D surface (set data2 = [Dx1 Dy1 Dz1 ; Dx2 Dy2 Dz2 ; ...])
% onto another 3D surface (set data1 = [dx1 dy1 dz1 ; dx2 dy2 dz2 ; ...])
%
% FORMAT [M1] = spm_eeg_inv_icp(data1,data2,fid1,fid2,Fmri,Fhsp,[aff])
% Input:
% data1      - locations of the first set of points corresponding to the
%            3D surface to register onto (p points)
% data2      - locations of the second set of points corresponding to the
%            second 3D surface to be registered (m points)
% fid1       - sMRI fiducials
% fid2       - sens fiducials
% Fmri       - graphics handle for sMRI points
% Fhsp       - graphics handle for headshape
% aff        - flag for 12 - parameter affine transform
%
% Output:
% M1         - 4 x 4 affine transformation matrix for sensor space
%==========================================================================
% Adapted from (http://www.csse.uwa.edu.au/~ajmal/icp.m) written by Ajmal Saeed Mian {ajmal@csse.uwa.edu.au}
% Computer Science, The University of Western Australia.

% Jeremie Mattout & Guillaume Flandin

% Landmarks (fiduciales) based registration
% Fiducial coordinates must be given in the same order in both files

% use figure and fiducials if specified
%--------------------------------------------------------------------------
try, fid1; catch, fid1 = []; end
try, fid2; catch, fid2 = []; end
try, aff;  catch, aff  = 0;  end


% initialise rotation and translation of sensor space
%--------------------------------------------------------------------------
M1    = speye(4,4);
tri   = delaunayn(data1');
for i = 1:16

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
        % 6-parmaeter affine (i.e. rigid body)
        %----------------------------------------------------------------------
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

end
return
%==========================================================================


%==========================================================================
function [M1] = spm_eeg_inv_rigidreg(data1, data2)
M       = spm_detrend(data1');
S       = spm_detrend(data2');
[U A V] = svd(S'*M);
R1      = V*U';
if det(R1) < 0
    B      = eye(3);
    B(3,3) = det(V*U');
    R1     = V*B*U';
end
t1 = mean(data1,2) - R1*mean(data2,2);
M1 = [R1 t1; 0 0 0 1];

return

%==========================================================================

function [h] = headcoordinates(nas, lpa, rpa)

% HEADCOORDINATES returns the homogenous coordinate transformation matrix
% that converts the specified fiducials in any coordinate system (e.g. MRI)
% into the rotated and translated headccordinate system.
%
% [h] = headcoordinates(nas, lpa, rpa)
%
% The headcoordinate system in CTF is defined as follows:
% the origin is exactly between lpa and rpa
% the X-axis goes towards nas
% the Y-axis goes approximately towards lpa, orthogonal to X and in the plane spanned by the fiducials
% the Z-axis goes approximately towards the vertex, orthogonal to X and Y
%
% See also WARPING, WARP3D

% Copyright (C) 2003 Robert Oostenveld
%
% $Log: headcoordinates.m,v $
% Revision 1.1  2004/09/27 16:00:04  roboos
% initial submission
%

% ensure that they are row vectors
lpa = lpa(:)';
rpa = rpa(:)';
nas = nas(:)';

% compute the origin and direction of the coordinate axes in MRI coordinates

% follow CTF convention
origin = [lpa+rpa]/2;
dirx = nas-origin;
dirx = dirx/norm(dirx);
dirz = cross(dirx,lpa-rpa);
dirz = dirz/norm(dirz);
diry = cross(dirz,dirx);

% compute the rotation matrix
rot = eye(4);
rot(1:3,1:3) = inv(eye(3) / [dirx; diry; dirz]);
% compute the translation matrix
tra = eye(4);
tra(1:4,4)   = [-origin(:); 1];
% compute the full homogenous transformation matrix from these two
h = rot * tra;
