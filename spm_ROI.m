function [xY, XYZmm, j] = spm_ROI(xY, XYZmm)
% Region of Interest specification
% FORMAT xY = spm_ROI(xY)
% xY     - VOI structure
%    xY.def      - VOI definition [sphere, box, mask, cluster, all]
%    xY.rej      - cell array of disabled VOI definition options
%    xY.xyz      - centre of VOI {mm}
%    xY.spec     - VOI definition parameters
%    xY.str      - description of the VOI
%
% FORMAT [xY, XYZmm, j] = spm_ROI(xY, XYZmm)
% XYZmm  - [3xm] locations of voxels {mm}
%          If an image filename, an spm_vol structure or a NIfTI object is
%          given instead, XYZmm will be initialised to all voxels within
%          the field of view of that image.
%
% XYZmm  - [3xn] filtered locations of voxels {mm} (m>=n) within VOI xY
% j      - [1xn] indices of input locations XYZmm within VOI xY
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston, Guillaume Flandin
% $Id: spm_ROI.m 3960 2010-06-30 17:41:24Z ged $

if nargin < 2 && nargout > 1
    error('Too many output arguments.');
end

try, xY; catch, xY = []; end

%-Specify ROI
%==========================================================================
if ~isfield(xY,'def')
    def        = {'sphere','box','cluster','mask'};
    if isfield(xY,'rej')
        if ~isfield(xY,'M')
            xY.rej = {xY.rej{:} 'cluster'};
        end
    else
        if isfield(xY,'M')
            xY.rej = {};
        else
            xY.rej = {'cluster'};
        end
    end
    [q, i] = setdiff(def,xY.rej);
    def    = def(sort(i));
    xY.def     = spm_input('VOI definition...','!+1','b',def,[],1);
end

%-ROI parameters
%--------------------------------------------------------------------------
switch lower(xY.def)

    case 'sphere'
    %----------------------------------------------------------------------
    if ~isfield(xY,'xyz') || isempty(xY.xyz)
        xY.xyz = spm_input('sphere centre [x y z] {mm}',...
            '!+0','r','0 0 0',3);
    end
    if ~isfield(xY,'spec')
        xY.spec = spm_input('sphere radius (mm)','!+0','r',0,1,[0,Inf]);
    end
    xY.str = sprintf('%0.1fmm sphere',xY.spec);

    case 'box'
    %----------------------------------------------------------------------
    if ~isfield(xY,'xyz') || isempty(xY.xyz)
        xY.xyz = spm_input('box centre [x y z] {mm}',...
            '!+0','r','0 0 0',3);
    end
    if ~isfield(xY,'spec')
        xY.spec = spm_input('box dimensions [x y z] {mm}',...
            '!+0','r','0 0 0',3);
    end
    if length(xY.spec) < 3
        xY.spec = xY.spec(1)*[1 1 1];
    end
    xY.str = sprintf('%0.1f x %0.1f x %0.1f mm box',xY.spec);
    
    case 'mask'
    %----------------------------------------------------------------------
    if ~isfield(xY,'spec')
        xY.spec = spm_vol(spm_select(1,'image','Specify Mask'));
    else
        if ~isstruct(xY.spec)
            xY.spec = spm_vol(xY.spec);
        end
    end
    str    = spm_str_manip(xY.spec.fname,'a30x');
    xY.str = sprintf('image mask: %s',str); 
        
    case 'cluster'
    %----------------------------------------------------------------------
    if ~isfield(xY,'xyz') || isempty(xY.xyz)
        xY.xyz = spm_input('seed voxel [x y z] {mm}',...
            '!+0','r','0 0 0',3);
    end
    if ~isfield(xY,'M')
        xY.M = spm_input('affine transformation matrix',...
            '!+0','r','0 0 0',[4 4]);
    end
    xY.spec = [];
    xY.str  = sprintf('cluster (seed voxel: %0.1f %0.1f %0.1f)',xY.xyz);
    
    case 'all'
    %----------------------------------------------------------------------
    xY.str  = 'all';
    
    otherwise
    %----------------------------------------------------------------------
    error('Unknown VOI type.');
    
end

if nargin < 2, return; end

%-'Estimate' ROI
%==========================================================================

%-Argument check
%--------------------------------------------------------------------------
if ischar(XYZmm) && isempty(XYZmm)
    XYZmm = spm_select(1,'image','Specify Image');
end
if ischar(XYZmm), XYZmm = spm_vol(XYZmm); end
if isa(XYZmm,'nifti')
    XYZmm    = struct('dim',size(XYZmm.dat), 'mat',XYZmm.mat);
end
if isstruct(XYZmm) % spm_vol
    [R,C,P]  = ndgrid(1:XYZmm.dim(1),1:XYZmm.dim(2),1:XYZmm.dim(3));
    RCP      = [R(:)';C(:)';P(:)'];
    clear R C P
    RCP(4,:) = 1;
    XYZmm    = XYZmm.mat(1:3,:)*RCP;    
end
if isempty(XYZmm), XYZmm = zeros(3,0); end

%-Filter location of voxels
%--------------------------------------------------------------------------
Q          = ones(1,size(XYZmm,2));

switch lower(xY.def)

    case 'sphere'
    %----------------------------------------------------------------------
    j      = find(sum((XYZmm - xY.xyz*Q).^2) <= xY.spec^2);

    case 'box'
    %----------------------------------------------------------------------
    j      = find(all(abs(XYZmm - xY.xyz*Q) <= xY.spec(:)*Q/2));
    
    case 'mask'
    %----------------------------------------------------------------------
    XYZ    = xY.spec.mat \ [XYZmm; Q];
    j      = find(spm_sample_vol(xY.spec, XYZ(1,:), XYZ(2,:), XYZ(3,:),0) > 0);
    
    case 'cluster'
    %----------------------------------------------------------------------
    [x i]  = spm_XYZreg('NearestXYZ',xY.xyz,XYZmm);
    XYZ    = round(xY.M \ [XYZmm; Q]);
    A      = spm_clusters(XYZ);
    j      = find(A == A(i));
    
    case 'all'
    %----------------------------------------------------------------------
    j      = 1:size(XYZmm,2);
    
    otherwise
    %----------------------------------------------------------------------
    error('Unknown VOI type.');
    
end

XYZmm      = XYZmm(:,j);
if strcmpi(xY.def,'mask') && ~isempty(XYZmm), xY.xyz = mean(XYZmm,2); end
