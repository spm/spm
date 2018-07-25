function varargout = spm_surf(P,mode,thresh)
% Surface extraction
% FORMAT spm_surf(P,mode,thresh)
%
% P      - char array of filenames
%          Usually, this will be c1xxx.<ext> & c2xxx.<ext> - grey and white
%          matter segments created using the segmentation routine.
% mode   - operation mode [1: rendering, 2: surface, 3: both]
% thresh - vector or threshold values for extraction [default: 0.5]
%          This is only relevant for extracting surfaces, not rendering.
%
% Generated files (depending on 'mode'):
% A "render_xxx.mat" file can be produced that can be used for
% rendering activations on to, see spm_render.
%
% A "xxx.surf.gii" file can also be written, which is created using
% Matlab's isosurface function.
% This extracted brain surface can be viewed using code something like:
%    FV = gifti(spm_select(1,'mesh','Select surface data'));
%    FV = export(FV,'patch');
%    fg = spm_figure('GetWin','Graphics');
%    ax = axes('Parent',fg);
%    p  = patch(FV, 'Parent',ax,...
%           'FaceColor', [0.8 0.7 0.7], 'FaceVertexCData', [],...
%           'EdgeColor', 'none',...
%           'FaceLighting', 'gouraud',...
%           'SpecularStrength' ,0.7, 'AmbientStrength', 0.1,...
%           'DiffuseStrength', 0.7, 'SpecularExponent', 10);
%    set(0,'CurrentFigure',fg);
%    set(fg,'CurrentAxes',ax);
%    l  = camlight(-40, 20); 
%    axis image;
%    rotate3d on;
%
% FORMAT out = spm_surf(job)
% 
% Input
% A job structure with fields
% .data   - cell array of filenames
% .mode   - operation mode
% .thresh - thresholds for extraction
% Output
% A struct with fields (depending on operation mode)
% .rendfile - cellstring containing render filename
% .surffile - cellstring containing surface filename(s)
%__________________________________________________________________________
%
% This surface extraction is not particularly sophisticated.  It simply
% smooths the data slightly and extracts the surface at a threshold of
% 0.5. The input segmentation images can be manually cleaned up first using
% e.g., MRIcron.
%__________________________________________________________________________
% Copyright (C) 2002-2018 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_surf.m 7381 2018-07-25 10:27:54Z guillaume $

SVNrev = '$Rev: 7381 $';

spm('FnBanner',mfilename,SVNrev);
spm('FigName','Surface');

%-Get input: filenames 'P'
%--------------------------------------------------------------------------
try
    if isstruct(P)
        if all(isfield(P,{'data','mode','thresh'}))
            job    = P;
            P      = char(job.data);
            mode   = job.mode;
            thresh = job.thresh;
        else
            % P is an spm_vol struct array
        end
    end
catch
    [P, sts] = spm_select([1 Inf],'image','Select images');
    if ~sts, varargout = {}; return; end
end

%-Get input: operation mode 'mode'
%--------------------------------------------------------------------------
try
    mode;
catch
    mode = spm_input('Save','+1','m',...
        ['Save Rendering|'...
         'Save Extracted Surface|'...
         'Save Rendering and Surface'],[1 2 3],3);
end

%-Get input: threshold for extraction 'thresh'
%--------------------------------------------------------------------------
try
    thresh;
catch
    thresh = 0.5;
end

%-Surface extraction
%--------------------------------------------------------------------------
spm('FigName','Surface: working');
spm('Pointer','Watch');
out = surfext(P,mode,thresh);
spm('Pointer','Arrow');
spm('FigName','Surface: done');

if nargout > 0
    varargout{1} = out;
end


%==========================================================================
% function out = surfext(P,mode,thresh)
%==========================================================================
function out = surfext(P,mode,thresh)

V  = spm_vol(P);
br = zeros(V(1).dim(1:3));
for i=1:V(1).dim(3)
    B         = spm_matrix([0 0 i]);
    tmp       = spm_slice_vol(V(1),B,V(1).dim(1:2),1);
    for j=2:numel(V)
        M   = V(j).mat\V(1).mat*B;
        tmp = tmp + spm_slice_vol(V(j),M,V(1).dim(1:2),1);
    end
    br(:,:,i) = tmp;
end

% Build a 3x3x3 seperable smoothing kernel and smooth
%--------------------------------------------------------------------------
kx = [0.75 1 0.75];
ky = [0.75 1 0.75];
kz = [0.75 1 0.75];
sm = sum(kron(kron(kz,ky),kx))^(1/3);
kx = kx/sm; ky = ky/sm; kz = kz/sm;
spm_conv_vol(br,br,kx,ky,kz,-[1 1 1]);

[pth,nam] = spm_fileparts(V(1).fname);

if any(mode==[1 3])
    % Produce rendering
    %----------------------------------------------------------------------
    out.rendfile{1} = fullfile(pth,['render_' nam '.mat']);
    tmp = struct('dat',br,'dim',size(br),'mat',V(1).mat);
    renviews(tmp,out.rendfile{1});
end

if any(mode==[2 3])
    % Produce extracted surface
    %----------------------------------------------------------------------
    for k=1:numel(thresh)
        [faces,vertices] = isosurface(br,thresh(k));

        % Swap around x and y because isosurface does for some
        % weird and wonderful reason.
        Mat      = V(1).mat(1:3,:)*[0 1 0 0;1 0 0 0;0 0 1 0; 0 0 0 1];
        vertices = (Mat*[vertices' ; ones(1,size(vertices,1))])';
        if numel(thresh)==1
            nam1 = nam;
        else
            nam1 = sprintf('%s-%d',nam,k);
        end
        out.surffile{k} = fullfile(pth,[nam1 '.surf.gii']);
        save(gifti(struct('faces',faces,'vertices',vertices)),out.surffile{k});
    end
end


%==========================================================================
% function renviews(V,oname)
%==========================================================================
function renviews(V,oname)
% Produce images for rendering activations to
%
% FORMAT renviews(V,oname)
% V     - mapped image to render, or alternatively
%         a structure of:
%         V.dat - 3D array
%         V.dim - size of 3D array
%         V.mat - affine mapping from voxels to millimeters
% oname - the name of the render.mat file.
%__________________________________________________________________________
%
% Produces a matrix file "render_xxx.mat" which contains everything that
% "spm_render" is likely to need.
%
% Ideally, the input image should contain values in the range of zero
% and one, and be smoothed slightly.  A threshold of 0.5 is used to
% distinguish brain from non-brain.
%__________________________________________________________________________

fprintf('%-40s: ','Rendering...');                                      %-#
                                
rend{1} = make_struct(V,[pi 0 pi/2]);    % Transverse 1
rend{2} = make_struct(V,[0 0 pi/2]);     % Transverse 2
rend{3} = make_struct(V,[0 pi/2 pi]);    % Sagittal 1
rend{4} = make_struct(V,[0 pi/2 0]);     % Sagittal 2
rend{5} = make_struct(V,[pi/2 pi/2 0]);  % Coronal 1
rend{6} = make_struct(V,[pi/2 pi/2 pi]); % Coronal 2

fprintf('%30s\n','...done');                                            %-#

fprintf('Writing Render file in:\n  %s\n',spm_file(oname,'cpath'));     %-#
save(oname,'rend', spm_get_defaults('mat.format'));

if ~spm('CmdLine')
    disp_renderings(rend);
    spm_print;
end


%==========================================================================
% function str = make_struct(V,thetas)
%==========================================================================
function str = make_struct(V,thetas)
[D,M]     = matdim(V.dim(1:3),V.mat,thetas);
[ren,dep] = make_pic(V,M*V.mat,D);
str       = struct('M',M,'ren',ren,'dep',dep);


%==========================================================================
% function [ren,zbuf] = make_pic(V,M,D)
%==========================================================================
function [ren,zbuf] = make_pic(V,M,D)
% A bit of a hack to try and make spm_render_vol produce some slightly
% prettier output.  It kind of works...
if isfield(V,'dat'), vv = V.dat; else vv = V; end;
[REN, zbuf, X, Y, Z] = spm_render_vol(vv, M, D, [0.5 1]);
fw        = max(sqrt(sum(M(1:3,1:3).^2)));
msk       = find(zbuf==1024);
brn       = ones(size(X));
brn(msk)  = 0;
brn       = spm_conv(brn,fw);
X(msk)    = 0;
Y(msk)    = 0;
Z(msk)    = 0;
msk       = find(brn<0.5);
tmp       = brn;
tmp(msk)  = 100000;
sX        = spm_conv(X,fw)./tmp;
sY        = spm_conv(Y,fw)./tmp;
sZ        = spm_conv(Z,fw)./tmp;
zbuf      = spm_conv(zbuf,fw)./tmp;
zbuf(msk) = 1024;

vec       = [-1 1 3]; % The direction of the lighting.
vec       = vec/norm(vec);
[t,dx,dy,dz] = spm_sample_vol(vv,sX,sY,sZ,3);
IM        = inv(diag([0.5 0.5 1])*M(1:3,1:3))';
ren       = IM(1:3,1:3)*[dx(:)' ; dy(:)' ; dz(:)'];
len       = sqrt(sum(ren.^2,1))+eps;
ren       = [ren(1,:)./len ; ren(2,:)./len ; ren(3,:)./len];
ren       = reshape(vec*ren,[size(dx) 1]);
ren(ren<0) = 0;
ren(msk)  = ren(msk)-0.2;
ren       = ren*0.8+0.2;
mx        = max(ren(:));
ren       = ren/mx;


%==========================================================================
% function disp_renderings(rend)
%==========================================================================
function disp_renderings(rend)
Fgraph = spm_figure('GetWin','Graphics');
spm_results_ui('Clear',Fgraph);
hght = 0.95;
nrow = ceil(length(rend)/2);
ax   = axes('Parent',Fgraph,'units','normalized',...
    'Position',[0, 0, 1, hght],'Visible','off');
image(0,'Parent',ax);
set(ax,'YTick',[],'XTick',[]);

for i=1:length(rend)
    ren = rend{i}.ren;
    ax  = axes('Parent',Fgraph,'units','normalized',...
            'Position',[rem(i-1,2)*0.5, floor((i-1)/2)*hght/nrow, 0.5, hght/nrow],...
            'Visible','off');
    image(ren*64,'Parent',ax);
    set(ax,'DataAspectRatio',[1 1 1], ...
        'PlotBoxAspectRatioMode','auto',...
        'YTick',[],'XTick',[],'XDir','normal','YDir','normal');
end


%==========================================================================
% function [d,M] = matdim(dim,mat,thetas)
%==========================================================================
function [d,M] = matdim(dim,mat,thetas)
R = spm_matrix([0 0 0 thetas]);
bb = [[1 1 1];dim(1:3)];
c  = [  bb(1,1) bb(1,2) bb(1,3) 1
        bb(1,1) bb(1,2) bb(2,3) 1
        bb(1,1) bb(2,2) bb(1,3) 1
        bb(1,1) bb(2,2) bb(2,3) 1
        bb(2,1) bb(1,2) bb(1,3) 1
        bb(2,1) bb(1,2) bb(2,3) 1
        bb(2,1) bb(2,2) bb(1,3) 1
        bb(2,1) bb(2,2) bb(2,3) 1]';
tc = diag([2 2 1 1])*R*mat*c;
tc = tc(1:3,:)';
mx = max(tc);
mn = min(tc);
M  = spm_matrix(-mn(1:2))*diag([2 2 1 1])*R;
d  = ceil(abs(mx(1:2)-mn(1:2)))+1;
