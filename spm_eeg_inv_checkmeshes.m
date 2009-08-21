function [h_ctx,h_skl,h_slp] = spm_eeg_inv_checkmeshes(varargin)
% Generate the tesselated surfaces of the inner-skull and scalp from binary volumes.
%
% FORMAT [h_ctx,h_skl,h_slp] = spm_eeg_inv_checkmeshes(S)
% Input:
% S         - input data struct (optional)
% Output:
% h_ctx     - handle to cortex patch
% h_skl     - handle to skull patch
% h_slp     - handle to scalp patch
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Jeremie Mattout
% $Id: spm_eeg_inv_checkmeshes.m 3328 2009-08-21 17:00:32Z vladimir $


% initialise
%--------------------------------------------------------------------------
[D,val] = spm_eeg_inv_check(varargin{:});
try
    disp(D.inv{val}.mesh);
    mesh = spm_eeg_inv_transform_mesh(eye(4), D.inv{val}.mesh);
    Mctx   = mesh.tess_ctx;
    Miskl  = mesh.tess_iskull;
    Moskl  = mesh.tess_oskull;
    Mslp   = mesh.tess_scalp;
catch
    warndlg('please create meshes')
    return
end

% SPM graphics figure
%--------------------------------------------------------------------------
Fgraph  = spm_figure('GetWin','Graphics'); spm_figure('Clear',Fgraph)

% Cortex Mesh Display
%--------------------------------------------------------------------------
h_ctx   = patch('vertices',Mctx.vert,'faces',Mctx.face,'EdgeColor','b','FaceColor','b');
hold on

% Inner-skull Mesh Display
%--------------------------------------------------------------------------
h_iskl   = patch('vertices',Miskl.vert,'faces',Miskl.face,'EdgeColor','r','FaceColor','none');

% Outer-skull Mesh Display
%--------------------------------------------------------------------------
h_oskl   = patch('vertices',Moskl.vert,'faces',Moskl.face,'EdgeColor',[1 .5 .35],'FaceColor','none');

% Inner Scalp Mesh Display
%--------------------------------------------------------------------------
h_slp   = patch('vertices',Mslp.vert,'faces',Mslp.face,'EdgeColor',[1 .7 .55],'FaceColor','none');


pls=0.05:0.2:0.9;
N = nifti(mesh.sMRI);
d=size(N.dat);
pls = round(pls.*d(3));
for i=1:numel(pls),
    [x,y,z]=ndgrid(1:d(1),1:d(2),pls(i));
    f1 = N.dat(:,:,pls(i));
    M = N.mat;
    x1 = M(1,1)*x+M(1,2)*y+M(1,3)*z+M(1,4);
    y1 = M(2,1)*x+M(2,2)*y+M(2,3)*z+M(2,4);
    z1 = M(3,1)*x+M(3,2)*y+M(3,3)*z+M(3,4);

    s=surf(x1,y1,z1,f1);
    set(s,'EdgeColor','none')
end

pnt = mesh.fid.fid.pnt;
h_fid = plot3(pnt(:,1), pnt(:,2), pnt(:,3), '.c', 'MarkerSize', 30);
h_fid_txt = text(pnt(:,1), pnt(:,2), pnt(:,3), mesh.fid.fid.label);

axis image off;
colormap('gray');
view(-135,45);
rotate3d on
drawnow
hold off

