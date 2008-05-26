function M1 = spm_eeg_inv_datareg(S)
% Rigid registration of the MEEG data and sMRI spaces
% rigid co-registration
%           1: fiducials based (3 landmarks: nasion, left ear, right ear)
%           2: surface matching between sensor mesh and headshape
%           (starts with a type 1 registration)
%
% FORMAT M1 = spm_eeg_inv_datareg(S)
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
% S.template  - 1 - input is a template (for EEG)
%               0 - input is an individual head model
%               2 - input is a template (for MEG) - enforce uniform scaling
%
% S.useheadshape - 1 use headshape matching 0 - don't
%
%
% Output:
% M1 = homogenous transformation matrix
%
% If a template is used, the senor locations are transformed using an
% affine (rigid body) mapping.  If headshape locations are supplied
% this is generalized to a full twelve parameter affine mapping (n.b.
% this might not be appropriate for MEG data).
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id: spm_eeg_inv_datareg.m 1726 2008-05-26 16:45:55Z vladimir $


if nargin == 0 || ~isstruct(S)
    error('Input struct is required');
end

if ~isfield(S, 'sens')
    error('MEEG sensors are missing');
else
    sens = forwinv_convert_units(S.sens, 'mm');
end

if ~isfield(S, 'mrifid')
    error('MRI fiducials are missing');
else
    mrifid = forwinv_convert_units(S.mrifid, 'mm');
    nfid = size(mrifid.fid.pnt, 1);
end

if ~isfield(S, 'meegfid')
    error('MEEG fiducials are missing');
else
    meegfid = forwinv_convert_units(S.meegfid, 'mm');
    meegfid.fid.pnt = meegfid.fid.pnt(1:nfid, :);
    meegfid.fid.label = meegfid.fid.label(1:nfid);
end

if isfield(S, 'vol');
    ishead = 1;
    vol = forwinv_convert_units(S.vol, 'mm');
end

if ~isfield(S, 'template')
    S.template = 0;
end


% Estimate-apply rigid body transform to sensor space
%--------------------------------------------------------------------------
M1 = spm_eeg_inv_rigidreg(mrifid.fid.pnt', meegfid.fid.pnt');

meegfid = forwinv_transform_headshape(M1, meegfid);

if S.template

    % constatined affine transform
    %--------------------------------------------------------------------------
    aff   = S.template;
    for i = 1:16

        % scale
        %----------------------------------------------------------------------
        M       = pinv(meegfid.fid.pnt(:))*mrifid.fid.pnt(:);
        M       = sparse(1:4,1:4,[M M M 1]);

        meegfid = forwinv_transform_headshape(M, meegfid);

        M1      = M*M1;

        % and move
        %----------------------------------------------------------------------
        M       = spm_eeg_inv_rigidreg(mrifid.fid.pnt', meegfid.fid.pnt');

        meegfid = forwinv_transform_headshape(M, meegfid);

        M1      = M*M1;

    end
else
    aff = 0;
end


% Surface matching between the scalp vertices in MRI space and
% the headshape positions in data space
%--------------------------------------------------------------------------
if ishead && ~isempty(meegfid.pnt) && S.useheadshape

    headshape = meegfid.pnt;
    scalpvert = mrifid.pnt;

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
    meegfid = forwinv_transform_headshape(M, meegfid);
    M1        = M*M1;
end
