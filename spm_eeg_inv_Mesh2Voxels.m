function [D] = spm_eeg_inv_Mesh2Voxels(varargin)
 
%======================================================================
% Non-linear interpolation of a Mesh contrast into MNI Voxel space
% This routine is used to produce a 3D image canonical sMRI
% space (in voxel coordinates) from a cortical mesh (3D surface).
% This yields a .nifti image of the summary statistics of the cortical
% activity for the effect of interest. This image can then enter the
% classical SPM routines for statistical testing.
% The [non-negative] root mean square contrast is smoothed both on the mesh
% (using a graph Laplacian and then in voxel-space using a conventional
% Gaussian filter.
%
% FORMAT D = spm_eeg_inv_Mesh2Voxels(D,val)
% Input:
% D        - input data struct (optional)
%     D.inv{val}.contrast.smooth  = smoothing in mm [12]
%     D.inv{val}.contrast.display = display (spm_image) flag [0]
% Output:
% D        - data struct including the new files and parameters
%
%==========================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience
 
% Jeremie Mattout, Stefan Kiebel & Karl Friston
% $Id: spm_eeg_inv_Mesh2Voxels.m 459 2006-02-24 11:55:50Z jeremie $
 
% checks
%--------------------------------------------------------------------------
[D,val]  = spm_eeg_inv_check(varargin{:});
 
% display flag
%--------------------------------------------------------------------------
try
    Disp = D.inv{val}.contrast.display;
catch
    Disp = 0;
end
 
% Scaling factor for images
%--------------------------------------------------------------------------
try
    scalefactor = D.inv{val}.contrast.scalefactor;
catch
    scalefactor = 1e6;
end
 
% smoothing FWHM (mm)
%--------------------------------------------------------------------------
fprintf('writing and smoothing image - please wait\n')
try
    smoothparam = D.inv{val}.contrast.smooth;
catch
    smoothparam = 8;
    D.inv{val}.contrast.smooth = smoothparam;
end
 
 
% extract variables
%----------------------------------------------------------------------
sMRIfile = fullfile(spm('dir'),'templates','T2.nii');
vert     = D.inv{val}.mesh.tess_mni.vert;
face     = D.inv{val}.mesh.tess_mni.face;
 
% Get mesh
%--------------------------------------------------------------------------
Vin   = spm_vol(sMRIfile);
nd    = D.inv{val}.inverse.Nd;
nv    = size(vert,1);
nf    = size(face,1);
 
% Compute a densely sampled triangular mask
%--------------------------------------------------------------------------
[tx, ty] = meshgrid(0.05:0.1:0.95, 0.05:0.1:0.95);
tx    = tx(:);
ty    = ty(:);
ti    = find(sum([tx ty]') <= 0.9);
t     = [tx(ti) ty(ti)];
np    = length(t);
 
% Map the template (square) triangle onto each face of the Cortical Mesh
%--------------------------------------------------------------------------
P1    = vert(face(:,1),:);
P2    = vert(face(:,2),:);
P3    = vert(face(:,3),:);
 
If    = speye(nf);
alpha = t(:,1);
beta  = t(:,2);
teta  = ones(np,1) - alpha - beta;
clear t tx ty ti
 
DenseCortex = kron(If,alpha)*P2 + kron(If,beta)*P3 + kron(If,teta)*P1;
clear P1 P2 P3
 
% Transform the sampling point coordinates from mm to voxels
%--------------------------------------------------------------------------
VoxelCoord  = round( spm_eeg_inv_mm2vx(DenseCortex,Vin.mat) )';
clear DenseCortex
 
% Get Garph Laplacian for smoothing on the cortical surface
%--------------------------------------------------------------------------
A     = spm_eeg_inv_meshdist(vert,face,0);
GL    = speye(nd,nd) + (A - spdiags(sum(A,2),0,nd,nd))/16;

 
% Interpolate the values in each vertex to compute the values at each
% sampling point of the triangles (cycle over conditions)
%==========================================================================
[PTH,NAME,EXT] = fileparts(D.fname);
for c = 1:length(D.inv{val}.contrast.GW)
    
    
    % Smooth on the cortical surface
    %----------------------------------------------------------------------
    Contrast = D.inv{val}.contrast.GW{c};
    Contrast = sparse(D.inv{val}.inverse.Is,1,Contrast,nd,1);
    for i = 1:32
        Contrast = GL*Contrast;
    end
    
    Outputfilename = fullfile(D.path,sprintf(  'w_%s_%.0f_%.0f.nii',NAME,val,c));
    Outputsmoothed = fullfile(D.path,sprintf( 'sw_%s_%.0f_%.0f.nii',NAME,val,c));
    Outputlog      = fullfile(D.path,sprintf('lsw_%s_%.0f_%.0f.nii',NAME,val,c));
    InterpOp       = [teta alpha beta];
    SPvalues       = zeros(nf*np,1);
    Vout           = Vin;
    Vout.fname     = Outputfilename;
    Vout.dt(1)     = spm_type('int16');
    RECimage       = zeros(Vout.dim);
 
    % And interpolate those values into voxels
    %----------------------------------------------------------------------
    for i = 1:nf
        TextVal = Contrast(face(i,:));
        if any(TextVal)
            ValTemp        = InterpOp*TextVal;
            SPvalues( (i-1)*np+1 : i*np ) = sum(TextVal)*(ValTemp/sum(ValTemp));
            Vox            = VoxelCoord( (i-1)*np+1 : i*np , : );
            Val            = SPvalues( (i-1)*np+1 : i*np );
            [UnVox,I,J]    = unique(Vox,'rows');
            IndV           = sub2ind(Vin.dim,UnVox(:,1),UnVox(:,2),UnVox(:,3));
            K              = 1:length(J);
            MatV           = zeros(max(J),length(J));
            IndK           = sub2ind(size(MatV),J,K');
            MatV(IndK)     = Val;
            SV             = sum(MatV')';
            RECimage(IndV) = RECimage(IndV) + SV;
        end
    end
 
    % Write (scaled) image
    %----------------------------------------------------------------------
    RECimage  = RECimage*scalefactor;
    Vout      = spm_write_vol(Vout,RECimage);
    
    
    % Smoothing
    %----------------------------------------------------------------------
    spm_smooth(Vout,Outputsmoothed,smoothparam);
    str = 'Summary statisic image (RMS response) written:\n %s\n %s (smoothed)\n';
    fprintf(str,Outputfilename,Outputsmoothed)
    D.inv{val}.contrast.Vout{c}  = Vout;
    D.inv{val}.contrast.fname{c} = Outputsmoothed;   
 
end
 
% display
%==========================================================================
if Disp, spm_eeg_inv_image_display(D); end
