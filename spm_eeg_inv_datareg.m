function [varargout] = spm_eeg_inv_datareg(varargin)

%==========================================================================
% Rigid registration of the EEG/MEG data and sMRI spaces
%
% FORMAT D = spm_eeg_inv_datareg(S)
% rigid co-registration
%           1: fiducials based (3 landmarks: nasion, left ear, right ear)
%           2: surface matching between sensor mesh and headshape
%           (starts with a type 1 registration)
% Input:
% D             - input data struct (optional)
% Output:
% D             - data struct including the new files and parameters
%
% FORMAT [eeg2mri,sen_reg,fid_reg,hsp_reg,orient_reg,mri2eeg,hsp2eeg]
%          = spm_eeg_inv_datareg(sensors,fid_eeg,fid_mri,headshape,scalpvert,megorient,template)
% Input:
%
% sensors       - a matrix coordinate of the sensor
%               locations ([Sx1 Sy1 Sz1 ; Sx2 Sy2 Sz2 ; ...])
% fid_eeg       - the fiducial coordinates in sensor
%               space ([Nx Ny Nz ; LEx LEy LEz ; REx REy REz])
% fid_mri       - the fiducial coordinates in sMRI
%               space ([nx ny nz ; lex ley lez ; rex rey rez])
% headshape     - the headshape point coordinates in
%               sensor space ([hx1 hy1 hz1 ; hx2 hy2 hz2 ; ...])
% scalpvert     - the vertices coordinates of the scalp
%               tesselation in mri space ([vx1 vy1 vz1 ; vx2 vy2 vz2 ;
%               ...])
% megorient     - the sensor orientations ([ox1 oy1 oz1 ; ox2 oy2 oz2 ; ...])
%               (MEG only)
% template      - 0/1 switch to enable affine mri2eeg transform
%
% IMPORTANT: all the coordinates must be in the same units (mm).
%
% Output:
% eeg2mri       - rigid transformation (Rotation + Translation)
% sen_reg       - the registered sensor coordinates in sMRI space
% fid_reg       - the registered fiduicals coordinates in sMRI space
% hsp_reg       - the registered headshap point coordinates in sMRI space
% orient_reg    - the registrated sensor orientations (MEG only)
%
% If a template is used, the senor locations are transformed using an
% affine (rigid body) mapping.  If headshape locations are supplied 
% this is generalized to a full twelve parameter affine mapping (n.b.
% this might not be appropriate for MEG data).
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_datareg.m 734 2007-02-09 21:21:02Z karl $

% Check input arguments
%==========================================================================

% % spm_eeg_inv_datareg(D)
%--------------------------------------------------------------------------
if nargin < 3
    try
        [D val]   = spm_eeg_inv_check(varargin{:});
        sensors   = D.inv{val}.datareg.sensors;
        fid_eeg   = D.inv{val}.datareg.fid_eeg;
        fid_mri   = D.inv{val}.datareg.fid_mri;
        headshape = D.inv{val}.datareg.headshape;
        scalpvert = D.inv{val}.datareg.scalpvert;
        megorient = D.inv{val}.datareg.megorient;
        template  = D.inv{val}.mesh.template;
    catch
        D   = spm_eeg_inv_datareg_ui(varargin{:});
        D   = spm_eeg_inv_datareg(D);
    end

% spm_eeg_inv_datareg(sensors,fid_eeg,fid_mri,headshape,scalpvert,megorient,template)
%--------------------------------------------------------------------------
else
    sensors = varargin{1};
    fid_eeg = varargin{2};
    fid_mri = varargin{3};
    
    try, headshape = varargin{4}; catch, headshape = sparse(0,3); end
    try, scalpvert = varargin{5}; catch, scalpvert = sparse(0,3); end
    try, megorient = varargin{6}; catch, megorient = sparse(0,3); end
    try, template  = varargin{7}; catch, template  = 0;           end
end

% The fiducial coordinates must be in the same order (usually: NZ & LE, RE)
%--------------------------------------------------------------------------
nfid      = size(fid_eeg,1);
if  nfid ~= size(fid_mri,1)
    warndlg('Please specify the same number of MRI and EEG/MEG fiducials');
    return
end
    
% Estimate-apply rigid body transform to sensor space
%--------------------------------------------------------------------------
M1        = spm_eeg_inv_rigidreg(fid_mri', fid_eeg');
fid_eeg   = M1*[fid_eeg'; ones(1,nfid)];
fid_eeg   = fid_eeg(1:3,:)';


if template
    
    % constatined affine transform
    %--------------------------------------------------------------------------
    aff   = 1;
    for i = 1:16
        
        % scale
        %----------------------------------------------------------------------
        M       = pinv(fid_eeg(:))*fid_mri(:);
        M       = sparse(1:4,1:4,[M M M 1]);
        fid_eeg = M*[fid_eeg'; ones(1,nfid)];
        fid_eeg = fid_eeg(1:3,:)';
        M1      = M*M1;

        % and move
        %----------------------------------------------------------------------
        M       = spm_eeg_inv_rigidreg(fid_mri', fid_eeg');
        fid_eeg = M*[fid_eeg'; ones(1,nfid)];
        fid_eeg = fid_eeg(1:3,:)';
        M1      = M*M1;

    end
else    
    aff = 0;   
end

% assume headshape locations are registered to sensor fiducials
%--------------------------------------------------------------------------
M2   = M1;

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
    
    
    % move headshape locations (NB: future code will allow for hsp fiducials)
    %----------------------------------------------------------------------
    % fid_hsp = headshape(1:3,:);
    % M2      = spm_eeg_inv_rigidreg(fid_eeg', fid_hsp');
    headshape = M2*[headshape'; ones(1,size(headshape,1))];
    headshape = headshape(1:3,:)';

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
    headshape = M*[headshape'; ones(1,size(headshape,1))];
    headshape = headshape(1:3,:)';
    fid_eeg   = M*[fid_eeg'; ones(1,nfid)];
    fid_eeg   = fid_eeg(1:3,:)';
    M1        = M*M1;
    
end

% Update the sensor locations and orientation
%--------------------------------------------------------------------------
sensors   = M1*[sensors'; ones(1,size(sensors,1))];
sensors   = sensors(1:3,:)';
megorient = megorient*M1(1:3,1:3)';

% retain valid sensor locations for leadfield computation
%--------------------------------------------------------------------------
if nargin < 3
    try
        sens  = setdiff(D.channels.eeg, D.channels.Bad);
    catch
        sens  = D.channels.eeg;
        D.channels.Bad = [];
    end
    sensors   = sensors(sens,:);
end

% ensure sensors lie outside the scalp
%--------------------------------------------------------------------------
if length(scalpvert)
    for i = 1:4
        tri     = delaunayn(scalpvert);
        j       = dsearchn(scalpvert, tri, sensors);
        dist    = sqrt(sum(sensors.^2,2)./sum(scalpvert(j,:).^2,2));
        dist    = min(dist,1);
        sensors = diag(1./dist)*sensors;
    end
end

% Ouptut arguments
%--------------------------------------------------------------------------
if nargout == 1
    D.inv{val}.datareg.eeg2mri           = M1;
    D.inv{val}.datareg.sens_coreg        = sensors;
    D.inv{val}.datareg.fid_coreg         = fid_eeg;
    D.inv{val}.datareg.hsp_coreg         = headshape;
    D.inv{val}.datareg.sens_orient_coreg = megorient;
    D.inv{val}.datareg.sens              = sens;

    varargout{1} = D;
else
    
    % varargout = {RT,sensors_reg,fid_reg,headshape_reg,orient_reg}
    %----------------------------------------------------------------------
    varargout{1} = M1;
    varargout{2} = sensors;
    varargout{3} = fid_eeg;
    varargout{4} = headshape;
    varargout{5} = megorient;

end

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
