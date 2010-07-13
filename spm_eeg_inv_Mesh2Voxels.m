function [D] = spm_eeg_inv_Mesh2Voxels(varargin)
% converts a mesh-representation of M/EEG power into a smoothed image
% FORMAT [D] = spm_eeg_inv_Mesh2Voxels(D,val)
% Input:
% D        - input data struct (optional)
%
%     D.inv{val}.contrast.display = display (spm_image) flag [0]
%
% Output:
% D        - data struct including the new files and parameters and
%
%--------------------------------------------------------------------------
% Non-linear interpolation of a Mesh contrast into MNI Voxel space
% This routine is used to produce a 3D image canonical sMRI
% space (in voxel coordinates) from a cortical mesh (3D surface).
% This yields a .nifti image of the summary statistics of the cortical
% activity for the effect of interest. This image can then enter the
% classical SPM routines for statistical testing.
% The [non-negative] mean square contrast is smoothed both on the mesh
% (using a graph Laplacian and then in voxel-space using a conventional
% Gaussian filter.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_eeg_inv_Mesh2Voxels.m 3996 2010-07-13 22:20:40Z vladimir $

% checks
%--------------------------------------------------------------------------
[D,val]  = spm_eeg_inv_check(varargin{:});

% Switches
%--------------------------------------------------------------------------
try, Disp   = D.inv{val}.contrast.display;   catch, Disp   = 0;  end
try, space  = D.inv{val}.contrast.space;     catch, space  = 1;  end
try, smooth = D.inv{val}.contrast.smoothing; catch, smooth = 8; end


% smoothing FWHM (mm)
%--------------------------------------------------------------------------
fprintf('writing and smoothing image - please wait\n');                 %-#

% Get volume (MNI [1] or native [0] output image space)
%--------------------------------------------------------------------------
if space
    sMRIfile = fullfile(spm('dir'),'templates','T2.nii');
else
    sMRIfile = D.inv{val}.mesh.sMRI;
end
Vin   = spm_vol(sMRIfile);

% Tag (to identify particular contrast settings)
%--------------------------------------------------------------------------
woi   = D.inv{val}.contrast.woi;
foi   = D.inv{val}.contrast.fboi;
Nw    = size(woi,1);
for i = 1:Nw
    tag{i} = ['t' sprintf('%d_', woi(i,:)) 'f' sprintf('%d_', foi)];
end

% Get mesh
%--------------------------------------------------------------------------
if space
    m = D.inv{val}.mesh.tess_mni;
else
    m = export(gifti(D.inv{val}.mesh.tess_ctx),'spm');
end
vert  = m.vert;
face  = m.face;
nd    = D.inv{val}.inverse.Nd;
nv    = size(vert,1);
nf    = size(face,1);

% Compute a densely sampled triangular mask
%--------------------------------------------------------------------------
[tx ty] = meshgrid(0.05:0.1:0.95, 0.05:0.1:0.95);
tx    = tx(:);
ty    = ty(:);
ti    = find(sum([tx ty],2) <= 0.9);
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
VoxelCoord = Vin.mat\[DenseCortex';ones(1,size(DenseCortex,1))];
VoxelCoord = round(VoxelCoord(1:3,:)');
clear DenseCortex

% Get Graph Laplacian for smoothing on the cortical surface
%--------------------------------------------------------------------------
A     = spm_mesh_distmtx(struct('vertices',vert,'faces',face),0);
GL    = speye(nd,nd) + (A - spdiags(sum(A,2),0,nd,nd))/16;


% normalize and embed in 3-space
%==========================================================================
[PTH,NAME,EXT] = fileparts(D.fname);


% accumulate mean of log-contrasts (over trials)
%--------------------------------------------------------------------------
GW      = D.inv{val}.contrast.GW;

bytrial = iscell(GW{1});

Ne = [];
if bytrial
    for c = 1:numel(GW)
        Ne(c) = numel(GW{c});
    end
else
    Ne = ones(1, numel(GW));
end

Nj = numel(GW)/Nw;

k  = 1;
iw = [];
ie = [];
for c = 1:length(GW)
    if bytrial
        cGW = GW{c};
    else
        cGW = GW(c);
    end
    
    for t = 1:Ne(c)
        % Smooth on the cortical surface
        %----------------------------------------------------------------------
        ssq{k} = full(sparse(D.inv{val}.inverse.Is,1,cGW{t},nd,1));
        for i = 1:smooth
            ssq{k} = GL*ssq{k};
        end
        
        % compute (truncated) moment
        %----------------------------------------------------------------------
        lss        = log(ssq{k} + eps);
        i          = lss > (max(lss) - log(32));
        meanlss(k) = mean(lss(i));
        
        iw(k) = c;
        ie(k) = t;
        
        k = k+1;
    end
end

scale = exp(mean(meanlss));
for c = 1:numel(ssq)
    
    % Normalise
    %----------------------------------------------------------------------
    Contrast   = ssq{c}/scale;
    
    
    % Initialise image
    %----------------------------------------------------------------------
    con       = mod(iw(c) - 1, Nj) + 1;
    str       = tag{ceil(iw(c)/Nj)};    
    if bytrial
        fname     = fullfile(D.path,sprintf('%s_%.0f_%s%.0f_%.0f.nii',NAME,val,str,con, ie(c)));
    else
        fname     = fullfile(D.path,sprintf('%s_%.0f_%s%.0f.nii',NAME,val,str,con));
    end
    Vout           = struct(...
        'fname',     fname,...
        'dim',       Vin.dim,...
        'dt',        [spm_type('float32') spm_platform('bigend')],...
        'mat',       Vin.mat,...
        'pinfo',     [1 0 0]',...
        'descrip',   '');
    
    InterpOp  = [teta alpha beta];
    SPvalues  = zeros(nf*np,1);
    RECimage  = zeros(Vout.dim);
    
    
    % And interpolate those values into voxels
    %----------------------------------------------------------------------
    for i = 1:nf
        TextVal = Contrast(face(i,:));
        if any(TextVal)
            ValTemp        = InterpOp*TextVal;
            SPvalues((i-1)*np+1:i*np) = sum(TextVal)*(ValTemp/sum(ValTemp));
            Vox            = VoxelCoord( (i-1)*np+1 : i*np , : );
            Val            = SPvalues( (i-1)*np+1 : i*np );
            [UnVox,I,J]    = unique(Vox,'rows');
            IndV           = sub2ind(Vin.dim,UnVox(:,1),UnVox(:,2),UnVox(:,3));
            K              = 1:length(J);
            MatV           = zeros(max(J),length(J));
            IndK           = sub2ind(size(MatV),J,K');
            MatV(IndK)     = Val;
            SV             = sum(MatV,2);
            RECimage(IndV) = RECimage(IndV) + SV;
        end
    end
    
    % 3D smoothing and thresholding
    %----------------------------------------------------------------------
    spm_smooth(RECimage,RECimage,1);
    RECimage = RECimage.*(RECimage > exp(-8));
    
    % Write (smoothed and scaled) image
    %----------------------------------------------------------------------
    Vout  = spm_write_vol(Vout,RECimage);
    
    
    % save and report
    %----------------------------------------------------------------------
    str = 'Summary-statistic image written:\n %s\n';
    fprintf(str,fname);
    
    D.inv{val}.contrast.Vout{c}  = Vout;
    D.inv{val}.contrast.fname{c} = fname;
    
end

% display
%==========================================================================
if Disp, spm_eeg_inv_image_display(D); end
