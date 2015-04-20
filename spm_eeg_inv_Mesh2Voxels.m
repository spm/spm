function [D] = spm_eeg_inv_Mesh2Voxels(varargin)
% Convert a mesh representation of M/EEG power into a smoothed image
% FORMAT [D] = spm_eeg_inv_Mesh2Voxels(D,[val])
% Input:
% D        - MEEG object or filename of M/EEG mat-file (optional)
%
%     D.inv{val}.contrast.display:   display image at the end {true, [false]}
%     D.inv{val}.contrast.space:     native [0] or MNI {1} output space
%     D.inv{val}.contrast.format:    output file format {['image'], 'mesh'}
%     D.inv{val}.contrast.smoothing: # iterations for cortical smoothing
%
% Output:
% D        - MEEG object containing the new image filenames in fields:
%
%     D.inv{val}.contrast.fname
%__________________________________________________________________________
%
% Non-linear interpolation of a Mesh contrast into MNI Voxel space
% This routine is used to produce a 3D image canonical sMRI
% space (in voxel coordinates) from a cortical mesh (3D surface).
% This yields a NIfTI image of the summary statistics of the cortical
% activity for the effect of interest. This image can then enter the
% classical SPM routines for statistical testing.
% The [non-negative] mean square contrast is smoothed both on the mesh
% (using a graph Laplacian) and then in voxel-space using a conventional
% Gaussian filter.
%__________________________________________________________________________
% Copyright (C) 2007-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_eeg_inv_Mesh2Voxels.m 6412 2015-04-20 10:14:50Z vladimir $


SVNrev = '$Rev: 6412 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);

%-Parse input arguments
%--------------------------------------------------------------------------
[D,val]      = spm_eeg_inv_check(varargin{:});

%-Get options
%--------------------------------------------------------------------------
try, space   = D.inv{val}.contrast.space;     catch, space  = 1; end
try, fmt     = D.inv{val}.contrast.format;    catch, fmt    = 'image'; end
try, smooth  = D.inv{val}.contrast.smoothing; catch, smooth = 8; end
try, Disp    = D.inv{val}.contrast.display;   catch, Disp   = 0; end

%-Time and Frequency windows of interest
%--------------------------------------------------------------------------
woi          = D.inv{val}.contrast.woi;
foi          = D.inv{val}.contrast.fboi;
Nw           = size(woi,1);

%-Get output image field of view and resolution
%--------------------------------------------------------------------------
if space
    sMRIfile = fullfile(spm('dir'),'canonical','avg152T2.nii');
else
    sMRIfile = D.inv{val}.mesh.sMRI;
end
Vin          = spm_vol(sMRIfile);

%-Get surface mesh
%--------------------------------------------------------------------------
if space
    m        = export(gifti(D.inv{val}.mesh.tess_mni),'patch');
else
    m        = export(gifti(D.inv{val}.mesh.tess_ctx),'patch');
end
nd           = D.inv{val}.inverse.Nd;


%-Accumulate mean of log-contrasts (over trials)
%==========================================================================
GL      = spm_mesh_smooth(m);

GW      = D.inv{val}.contrast.GW;

bytrial = iscell(GW{1});

if bytrial
    Ne  = cellfun(@numel,GW);
else
    Ne  = ones(1, numel(GW));
end

Nj      = numel(GW)/Nw;

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
        %-Smooth on the cortical surface
        %------------------------------------------------------------------
        ssq{k} = full(sparse(D.inv{val}.inverse.Is,1,cGW{t},nd,1));
        ssq{k} = spm_mesh_smooth(GL,ssq{k},smooth);
        
        %-Compute (truncated) moment
        %------------------------------------------------------------------
        lss        = log(ssq{k} + eps);
        i          = lss > (max(lss) - log(32));
        meanlss(k) = mean(lss(i));
        
        iw(k) = c;
        ie(k) = t;
        
        k = k + 1;
    end
end

scale = exp(mean(meanlss));

%-Save as meshes
%==========================================================================
if strcmpi(fmt,'mesh')
    [pth,name] = fileparts(D.fname);
    tag = cell(1,Nw);
    for i = 1:Nw
        tag{i} = ['t' sprintf('%g_', woi(i,:)) 'f' sprintf('%d_', foi)];
    end
    
    %-Save mesh topology
    %----------------------------------------------------------------------
    g = gifti(m);
    save(g,fullfile(D.path,[name '.surf.gii']));
    
    %-Save mesh data
    %----------------------------------------------------------------------
    spm_progress_bar('Init',numel(ssq),'Exporting as meshes','');
    for c = 1:numel(ssq)
        
        fprintf('%s%30s',repmat(sprintf('\b'),1,30),...
            sprintf('...mesh %d/%d',c,numel(ssq)));                     %-#
        
        %-Filename
        %------------------------------------------------------------------
        con       = mod(iw(c) - 1, Nj) + 1;
        str       = tag{ceil(iw(c)/Nj)};
        if bytrial, bt = sprintf('_%.0f',ie(c)); else bt = ''; end
        fname     = fullfile(D.path,...
            sprintf('%s_%.0f_%s%.0f%s.gii', name, val, str, con, bt));
        
        %-Normalise
        %------------------------------------------------------------------
        Contrast = ssq{c} / scale;
        %Contrast = Contrast.*(Contrast > exp(-8));
        
        %-Write mesh
        %------------------------------------------------------------------
        g = gifti;
        g.cdata = Contrast;
        g.private.metadata(1).name = 'SurfaceID';
        g.private.metadata(1).value = [name '.surf.gii'];
        save(g, fname, 'ExternalFileBinary');
        
        %-Store filename
        %------------------------------------------------------------------
        D.inv{val}.contrast.fname{c} = fname;
        
        spm_progress_bar('Set', c);
        
    end
    
    spm_progress_bar('Clear');
    fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done');           %-#
    
    return;
end

%-Normalise and embed in 3D-space
%==========================================================================
fprintf('%-40s: %30s','Writing images','...please wait');               %-#

[pth,name] = fileparts(D.fname);
tag = cell(1,Nw);
for i = 1:Nw
    tag{i} = ['t' sprintf('%d_', woi(i,:)) 'f' sprintf('%d_', foi)];
end

spm_progress_bar('Init',numel(ssq),'Interpolating images','');

con       = mod(iw - 1, Nj) + 1;
win       = ceil(iw/Nj);

ucon      = unique(con);
uwin      = unique(win);

n = 0;
for c = 1:length(ucon)
    for w = 1:numel(uwin)
        
        ind = find((con == ucon(c)) & (win == uwin(w)));
        str = tag{uwin(w)};
        
        fname     = fullfile(D.path,...
            sprintf('%s_%.0f_%s%.0f.nii', name, val, str, ucon(c)));
        
        %-Initialise image
        %----------------------------------------------------------------------
        N      = nifti;
        N.dat  = file_array(fname, [Vin.dim, length(ind)], 'FLOAT32-LE');
        N.mat  = Vin.mat;
        N.mat0 = Vin.mat;
        create(N);
        
        
        for i = 1:length(ind)
            
            n = n+1;
            
            fprintf('%s%30s',repmat(sprintf('\b'),1,30),...
                sprintf('...image %d/%d',n,numel(ssq)));                        %-#            
            
            %-Normalise
            %----------------------------------------------------------------------
            Contrast = ssq{ind(i)} / scale;
            
            %-Interpolate those values into voxels
            %----------------------------------------------------------------------
            RECimage = spm_mesh_to_grid(m, Vin, Contrast);
            
            %-3D smoothing and thresholding
            %----------------------------------------------------------------------
            spm_smooth(RECimage, RECimage, 1);
            RECimage = RECimage.*(RECimage > exp(-8));
            
            %-Write (smoothed and scaled) image
            %----------------------------------------------------------------------            
            N.dat(:, :, :, i) = RECimage;
            
            %-Store filename
            %----------------------------------------------------------------------
            D.inv{val}.contrast.fname{n} = fname;
            
            spm_progress_bar('Set', n);
            
        end
    end
end

spm_progress_bar('Clear');
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done');               %-#

%-Display
%==========================================================================
if Disp && ~spm('CmdLine'), spm_eeg_inv_image_display(D); end
