function [M1, sens, meegfid] = spm_eeg_inv_datareg_eeg(S)
% Rigid registration of the EEG data and sMRI spaces
% rigid co-registration
%           1: fiducials based (3 landmarks: nasion, left ear, right ear)
%           2: surface matching between sensor mesh and headshape
%           (starts with a type 1 registration)
%
% FORMAT [M1, sens, fid] = spm_eeg_inv_datareg_eeg(S)
%
% Input:
%
% S  - input struct
% fields of S:
%
% S.sens - EEG sensors (struct)
% S.meegfid  - EEG fiducials (struct)
% S.vol - volume model
% S.mrifid = MRI fiducials
% S.template  - 1 - input is a template /0 - input is an individual head
%                   model
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
% $Id: spm_eeg_inv_datareg_eeg.m 1488 2008-04-27 14:11:48Z vladimir $


if nargin == 0 || ~isstruct(S)
    error('Input struct is required');
end

if ~isfield(S, 'sens')
    error('EEG sensors are missing');
else
    sens = convert_units(S.sens, 'mm');
end

if ~isfield(S, 'meegfid')
    error('EEG fiducials are missing');
else
    meegfid = convert_units(S.meegfid, 'mm');
end

if ~isfield(S, 'mrifid')
    error('MRI fiducials are missing');
else
    mrifid = convert_units(S.mrifid, 'mm');
end

if isfield(S, 'vol');
    ishead = 1;
    vol = convert_units(S.vol, 'mm');
end

if ~isfield(S, 'template')
    S.template = 0;
end


% Estimate-apply rigid body transform to sensor space
%--------------------------------------------------------------------------
M1 = spm_eeg_inv_rigidreg(mrifid.fid.pnt', meegfid.fid.pnt');

meegfid = transform_headshape(M1, meegfid);

if S.template

    % constatined affine transform
    %--------------------------------------------------------------------------
    aff   = 1;
    for i = 1:16

        % scale
        %----------------------------------------------------------------------
        M       = pinv(meegfid.fid.pnt(:))*mrifid.fid.pnt(:);
        M       = sparse(1:4,1:4,[M M M 1]);

        meegfid = transform_headshape(M, meegfid);

        M1      = M*M1;

        % and move
        %----------------------------------------------------------------------
        M       = spm_eeg_inv_rigidreg(mrifid.fid.pnt', meegfid.fid.pnt');

        meegfid = transform_headshape(M, meegfid);

        M1      = M*M1;

    end
else
    aff = 0;
end


% Surface matching between the scalp vertices in MRI space and
% the headshape positions in data space
%--------------------------------------------------------------------------
if ishead && ~isempty(meegfid.pnt)

    headshape = meegfid.pnt;
    scalpvert = vol.bnd(1).pnt;

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
    M    = spm_eeg_inv_icp(scalpvert',headshape',mrifid.fid.pnt',meegfid.fid.pnt',Fmri,Fhsp,aff);

    % transform headshape and eeg fiducials
    %----------------------------------------------------------------------
    meegfid = transform_headshape(M, meegfid);
    M1        = M*M1;
end


% Update the sensor locations and orientation
%--------------------------------------------------------------------------
sens = transform_sens(M1, sens);

return

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
