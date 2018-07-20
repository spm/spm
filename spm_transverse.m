function spm_transverse(varargin)
% Rendering of regional effects [SPM{T/F}] on transverse sections
% FORMAT spm_transverse('set',SPM,hReg)
% FORMAT spm_transverse('setcoords',xyzmm)
% FORMAT spm_transverse('clear')
%
% SPM    - structure containing SPM, distribution & filtering details
%          about the excursion set (xSPM)
%        - required fields are:
% .Z     - minimum of n Statistics {filtered on u and k}
% .STAT  - distribution {Z, T, X or F}     
% .u     - height threshold
% .XYZ   - location of voxels {voxel coords}
% .iM    - mm  -> voxels matrix
% .VOX   - voxel dimensions {mm}
% .DIM   - image dimensions {voxels}
%
% hReg   - handle of MIP XYZ registry object (see spm_XYZreg for details)
%
% spm_transverse automatically updates its co-ordinates from the
% registry, but clicking on the slices has no effect on the registry.
% i.e., the updating is one way only.
%
% See also: spm_getSPM
%__________________________________________________________________________
%
% spm_transverse is called by the SPM results section and uses
% variables in SPM and SPM to create three transverse sections though a
% background image.  Regional foci from the selected SPM{T/F} are
% rendered on this image.
%
% Although the SPM{.} adopts the neurological convention (left = left)
% the rendered images follow the same convention as the original data.
%__________________________________________________________________________
% Copyright (C) 1994-2018 Wellcome Trust Centre for Neuroimaging

% Karl Friston & John Ashburner
% $Id: spm_transverse.m 7376 2018-07-20 10:30:59Z guillaume $


switch lower(varargin{1})

    case 'set'
    % draw slices
    %----------------------------------------------------------------------
    init(varargin{2},varargin{3});

    case 'setcoords'
    % reposition
    %----------------------------------------------------------------------
    disp('Reposition');

    case 'clear'
    % clear
    %----------------------------------------------------------------------
    clear_global;
end

return;


%==========================================================================
% function init(SPM,hReg)
%==========================================================================
function init(SPM,hReg)

%-Get figure handles
%--------------------------------------------------------------------------
Fgraph = spm_figure('GetWin','Graphics');

%-Get the image on which to render
%--------------------------------------------------------------------------
[spms, sts]   = spm_select(1,'image','Select image for rendering on');
if ~sts, return; end
spm('Pointer','Watch');

%-Delete previous axis and their pagination controls (if any)
%--------------------------------------------------------------------------
spm_results_ui('Clear',Fgraph);

global transv
transv      = struct('blob',[],'V',spm_vol(spms),'h',[],'hReg',hReg,'fig',Fgraph);
transv.blob = struct('xyz', round(SPM.XYZ), 't',SPM.Z, 'dim',SPM.DIM(1:3),...
                     'iM',SPM.iM,...
                     'vox', sqrt(sum(SPM.M(1:3,1:3).^2)), 'u', SPM.u);

%-Get current location and convert to pixel co-ordinates
%--------------------------------------------------------------------------
xyzmm  = spm_XYZreg('GetCoords',transv.hReg);
xyz    = round(transv.blob.iM(1:3,:)*[xyzmm; 1]);
try
    units = SPM.units;
catch
    units = {'mm' 'mm' 'mm'};
end

%-Extract data from SPM [at one plane separation] and get background slices
%--------------------------------------------------------------------------
dim    = ceil(transv.blob.dim(1:3)'.*transv.blob.vox);
A      = transv.blob.iM*transv.V.mat;
hld    = 0;

zoomM  = inv(spm_matrix([0 0 -1  0 0 0  transv.blob.vox([1 2]) 1]));
zoomM1 =     spm_matrix([0 0  0  0 0 0  transv.blob.vox([1 2]) 1]);

Q      = find(abs(transv.blob.xyz(3,:) - xyz(3)) < 0.5);
T2     = full(sparse(transv.blob.xyz(1,Q),transv.blob.xyz(2,Q),transv.blob.t(Q),transv.blob.dim(1),transv.blob.dim(2)));
T2     = spm_slice_vol(T2,zoomM,dim([1 2]),[hld NaN]);
Q      = find(T2==0) ; T2(Q) = NaN;
D      = zoomM1*[1 0 0 0;0 1 0 0;0 0 1 -xyz(3);0 0 0 1]*A;
D2     = spm_slice_vol(transv.V,inv(D),dim([1 2]),1);
maxD   = max([max(D2(:)) eps]);
minD   = min([min(D2(:)) eps]);

if transv.blob.dim(3) > 1

    Q      = find(abs(transv.blob.xyz(3,:) - xyz(3)+1) < 0.5);
    T1     = full(sparse(transv.blob.xyz(1,Q),...
            transv.blob.xyz(2,Q),transv.blob.t(Q),transv.blob.dim(1),transv.blob.dim(2)));
    T1     = spm_slice_vol(T1,zoomM,dim([1 2]),[hld NaN]);
    Q      = find(T1==0) ; T1(Q) = NaN;
    D      = zoomM1*[1 0 0 0;0 1 0 0;0 0 1 -xyz(3)+1;0 0 0 1]*A;
    D1     = spm_slice_vol(transv.V,inv(D),dim([1 2]),1);
    maxD   = max([maxD ; D1(:)]);
    minD   = min([minD ; D1(:)]);

    Q      = find(abs(transv.blob.xyz(3,:) - xyz(3)-1) < 0.5);
    T3     = full(sparse(transv.blob.xyz(1,Q),...
            transv.blob.xyz(2,Q),transv.blob.t(Q),transv.blob.dim(1),transv.blob.dim(2)));
    T3     = spm_slice_vol(T3,zoomM,dim([1 2]),[hld NaN]);
    Q      = find(T3==0) ; T3(Q) = NaN;
    D      = zoomM1*[1 0 0 0;0 1 0 0;0 0 1 -xyz(3)-1;0 0 0 1]*A;
    D3     = spm_slice_vol(transv.V,inv(D),dim([1 2]),1);
    maxD   = max([maxD ; D3(:)]);
    minD   = min([minD ; D3(:)]);
end

mx     = max([max(T2(:)) eps]);
mn     = min([min(T2(:)) 0]);
D2     = (D2-minD)/(maxD-minD);
if transv.blob.dim(3) > 1,
    D1 = (D1-minD)/(maxD-minD);
    D3 = (D3-minD)/(maxD-minD);
    mx = max([mx ; T1(:) ; T3(:) ; eps]);
    mn = min([mn ; T1(:) ; T3(:) ; 0]);
end;

%-Configure {128 level} colormap
%--------------------------------------------------------------------------
cmap   = get(Fgraph,'Colormap');
if size(cmap,1) ~= 128
    figure(Fgraph);
    spm_colourmap('gray-hot');
    cmap = get(Fgraph,'Colormap');
end

D      = length(cmap)/2;
Q      = find(T2(:) > transv.blob.u); T2 = (T2(Q)-mn)/(mx-mn); D2(Q) = 1+1.51/D + T2; T2 = D*D2;

if transv.blob.dim(3) > 1
    Q  = find(T1(:) > transv.blob.u); T1 = (T1(Q)-mn)/(mx-mn); D1(Q) = 1+1.51/D + T1; T1 = D*D1;
    Q  = find(T3(:) > transv.blob.u); T3 = (T3(Q)-mn)/(mx-mn); D3(Q) = 1+1.51/D + T3; T3 = D*D3;
end

set(Fgraph,'Units','pixels')
siz    = get(Fgraph,'Position');
siz    = siz(3:4);

P = xyz.*transv.blob.vox';

%-Render activation foci on background images
%--------------------------------------------------------------------------
if transv.blob.dim(3) > 1
    zm     = min([(siz(1) - 120)/(dim(1)*3),(siz(2)/2 - 60)/dim(2)]);
    xo     = (siz(1)-(dim(1)*zm*3)-120)/2;
    yo     = (siz(2)/2 - dim(2)*zm - 60)/2;

    transv.h(1) = axes('Units','pixels','Parent',Fgraph,'Position',[20+xo 20+yo dim(1)*zm dim(2)*zm]);
    transv.h(2) = image(rot90(spm_grid(T1)),'Parent',transv.h(1));
    axis image; axis off;
    tmp = SPM.iM\[xyz(1:2)' (xyz(3)-1) 1]';
    
    ax=transv.h(1);tpoint=get(ax,'title');
    str=sprintf('z = %0.0f%s',tmp(3),units{3});
    set(tpoint,'string',str);
    
    transv.h(3) = line([1 1]*P(1),[0 dim(2)],'Color','w','Parent',transv.h(1));
    transv.h(4) = line([0 dim(1)],[1 1]*(dim(2)-P(2)+1),'Color','w','Parent',transv.h(1));

    transv.h(5) = axes('Units','pixels','Parent',Fgraph,'Position',[40+dim(1)*zm+xo 20+yo dim(1)*zm dim(2)*zm]);
    transv.h(6) = image(rot90(spm_grid(T2)),'Parent',transv.h(5));
    axis image; axis off;
    tmp = SPM.iM\[xyz(1:2)' xyz(3) 1]';
    
    ax=transv.h(5);tpoint=get(ax,'title');
    str=sprintf('z = %0.0f%s',tmp(3),units{3});
    set(tpoint,'string',str);
    
    transv.h(7) = line([1 1]*P(1),[0 dim(2)],'Color','w','Parent',transv.h(5));
    transv.h(8) = line([0 dim(1)],[1 1]*(dim(2)-P(2)+1),'Color','w','Parent',transv.h(5));

    transv.h(9) = axes('Units','pixels','Parent',Fgraph,'Position',[60+dim(1)*zm*2+xo 20+yo dim(1)*zm dim(2)*zm]);
    transv.h(10) = image(rot90(spm_grid(T3)),'Parent',transv.h(9));
    axis image; axis off;
    tmp = SPM.iM\[xyz(1:2)' (xyz(3)+1) 1]';

    ax=transv.h(9);tpoint=get(ax,'title');
    str=sprintf('z = %0.0f%s',tmp(3),units{3});
    set(tpoint,'string',str);
    
    transv.h(11) = line([1 1]*P(1),[0 dim(2)],'Color','w','Parent',transv.h(9));
    transv.h(12) = line([0 dim(1)],[1 1]*(dim(2)-P(2)+1),'Color','w','Parent',transv.h(9));
    
    % colorbar
    %----------------------------------------------------------------------
    q      = [80+dim(1)*zm*3+xo 20+yo 20 dim(2)*zm];
    if SPM.STAT=='P'
        str='Effect size';
    else
        str=[SPM.STAT ' value'];
    end
    transv.h(13) = axes('Units','pixels','Parent',Fgraph,'Position',q,'Visible','off');
    transv.h(14) = image([0 mx/32],[mn mx],(1:D)' + D,'Parent',transv.h(13));

    ax=transv.h(13);
    tpoint=get(ax,'title');
    set(tpoint,'string',str);
    set(tpoint,'FontSize',9);
    %title(ax,str,'FontSize',9);
    set(ax,'XTickLabel',[]);
    axis(ax,'xy');

else
    zm     = min([(siz(1) - 80)/dim(1),(siz(2)/2 - 60)/dim(2)]);
    xo     = (siz(1)-(dim(1)*zm)-80)/2;
    yo     = (siz(2)/2 - dim(2)*zm - 60)/2;

    transv.h(1) = axes('Units','pixels','Parent',Fgraph,'Position',[20+xo 20+yo dim(1)*zm dim(2)*zm]);
    transv.h(2) = image(rot90(spm_grid(T2)),'Parent',transv.h(1));
    axis image; axis off;
    title(sprintf('z = %0.0f%s',xyzmm(3),units{3}));
    transv.h(3) = line([1 1]*P(1),[0 dim(2)],'Color','w','Parent',transv.h(1));
    transv.h(4) = line([0 dim(1)],[1 1]*(dim(2)-P(2)+1),'Color','w','Parent',transv.h(1));
    
    % colorbar
    %----------------------------------------------------------------------
    q      = [40+dim(1)*zm+xo 20+yo 20 dim(2)*zm];
    transv.h(5) = axes('Units','pixels','Parent',Fgraph,'Position',q,'Visible','off');
    transv.h(6) = image([0 mx/32],[mn mx],(1:D)' + D,'Parent',transv.h(5));
    if SPM.STAT=='P'
        str='Effect size';
    else
        str=[SPM.STAT ' value'];
    end
    
    title(str,'FontSize',9);
    set(gca,'XTickLabel',[]);
    axis xy;
    
    
end

spm_XYZreg('Add2Reg',transv.hReg,transv.h(1), 'spm_transverse');

for h=transv.h,
    set(h,'DeleteFcn',@clear_global);
end

%-Reset pointer
%--------------------------------------------------------------------------
spm('Pointer','Arrow')
return;

%==========================================================================
% function reposition(xyzmm)
%==========================================================================
function reposition(xyzmm)
global transv
if ~isstruct(transv), return; end;

spm('Pointer','Watch');


%-Get current location and convert to pixel co-ordinates
%--------------------------------------------------------------------------
% xyzmm  = spm_XYZreg('GetCoords',transv.hReg)
xyz    = round(transv.blob.iM(1:3,:)*[xyzmm; 1]);

% extract data from SPM [at one plane separation]
% and get background slices
%--------------------------------------------------------------------------
dim    = ceil(transv.blob.dim(1:3)'.*transv.blob.vox);
A      = transv.blob.iM*transv.V.mat;
hld    = 0;

zoomM  = inv(spm_matrix([0 0 -1  0 0 0  transv.blob.vox([1 2]) 1]));
zoomM1 =     spm_matrix([0 0  0  0 0 0  transv.blob.vox([1 2]) 1]);

Q      = find(abs(transv.blob.xyz(3,:) - xyz(3)) < 0.5);
T2     = full(sparse(transv.blob.xyz(1,Q),transv.blob.xyz(2,Q),transv.blob.t(Q),transv.blob.dim(1),transv.blob.dim(2)));
T2     = spm_slice_vol(T2,zoomM,dim([1 2]),[hld NaN]);
Q      = find(T2==0) ; T2(Q) = NaN;
D      = zoomM1*[1 0 0 0;0 1 0 0;0 0 1 -xyz(3);0 0 0 1]*A;
D2     = spm_slice_vol(transv.V,inv(D),dim([1 2]),1);
maxD   = max([max(D2(:)) eps]);
minD   = min([min(D2(:)) 0]);

if transv.blob.dim(3) > 1

    Q      = find(abs(transv.blob.xyz(3,:) - xyz(3)+1) < 0.5);
    T1     = full(sparse(transv.blob.xyz(1,Q),...
            transv.blob.xyz(2,Q),transv.blob.t(Q),transv.blob.dim(1),transv.blob.dim(2)));
    T1     = spm_slice_vol(T1,zoomM,dim([1 2]),[hld NaN]);
    Q      = find(T1==0) ; T1(Q) = NaN;
    D      = zoomM1*[1 0 0 0;0 1 0 0;0 0 1 -xyz(3)+1;0 0 0 1]*A;
    D1     = spm_slice_vol(transv.V,inv(D),dim([1 2]),1);
    maxD   = max([maxD ; D1(:)]);
    minD   = min([minD ; D1(:)]);

    Q      = find(abs(transv.blob.xyz(3,:) - xyz(3)-1) < 0.5);
    T3     = full(sparse(transv.blob.xyz(1,Q),...
            transv.blob.xyz(2,Q),transv.blob.t(Q),transv.blob.dim(1),transv.blob.dim(2)));
    T3     = spm_slice_vol(T3,zoomM,dim([1 2]),[hld NaN]);
    Q      = find(T3==0) ; T3(Q) = NaN;
    D      = zoomM1*[1 0 0 0;0 1 0 0;0 0 1 -xyz(3)-1;0 0 0 1]*A;
    D3     = spm_slice_vol(transv.V,inv(D),dim([1 2]),1);
    maxD   = max([maxD ; D3(:)]);
    minD   = min([minD ; D3(:)]);
end

mx     = max([max(T2(:)) eps]);
mn     = min([min(T2(:)) 0]);
D2     = (D2-minD)/(maxD-minD);
if transv.blob.dim(3) > 1,
    D1 = (D1-minD)/(maxD-minD);
    D3 = (D3-minD)/(maxD-minD);
    mx = max([mx ; T1(:) ; T3(:) ; eps]);
    mn = min([mn ; T1(:) ; T3(:) ; 0]);
end;

%-Configure {128 level} colormap
%--------------------------------------------------------------------------
cmap   = get(transv.fig,'Colormap');
if size(cmap,1) ~= 128
    figure(transv.fig);
    spm_colourmap('gray-hot');
    cmap = get(transv.fig,'Colormap');
end

D      = length(cmap)/2;
Q      = find(T2(:) > transv.blob.u); T2 = (T2(Q)-mn)/(mx-mn); D2(Q) = 1+1.51/D + T2; T2 = D*D2;

if transv.blob.dim(3) > 1
    Q  = find(T1(:) > transv.blob.u); T1 = (T1(Q)-mn)/(mx-mn); D1(Q) = 1+1.51/D + T1; T1 = D*D1;
    Q  = find(T3(:) > transv.blob.u); T3 = (T3(Q)-mn)/(mx-mn); D3(Q) = 1+1.51/D + T3; T3 = D*D3;
end

P = xyz.*transv.blob.vox';

%-Render activation foci on background images
%--------------------------------------------------------------------------
if transv.blob.dim(3) > 1

    set(transv.h(2),'Cdata',rot90(spm_grid(T1)));
    tmp = transv.blob.iM\[xyz(1:2)' (xyz(3)-1) 1]';
    set(get(transv.h(1),'Title'),'String',sprintf('z = %0.0fmm',tmp(3)));
    set(transv.h(3),'Xdata',[1 1]*P(1),'Ydata',[0 dim(2)]);
    set(transv.h(4),'Xdata',[0 dim(1)],'Ydata',[1 1]*(dim(2)-P(2)+1));

    set(transv.h(6),'Cdata',rot90(spm_grid(T2)));
    set(get(transv.h(5),'Title'),'String',sprintf('z = %0.0fmm',xyzmm(3)));
    set(transv.h(7),'Xdata',[1 1]*P(1),'Ydata',[0 dim(2)]);
    set(transv.h(8),'Xdata',[0 dim(1)],'Ydata',[1 1]*(dim(2)-P(2)+1));

    set(transv.h(10),'Cdata',rot90(spm_grid(T3)));
    tmp = transv.blob.iM\[xyz(1:2)' (xyz(3)+1) 1]';
    set(get(transv.h(9),'Title'),'String',sprintf('z = %0.0fmm',tmp(3)));
    set(transv.h(11),'Xdata',[1 1]*P(1),'Ydata',[0 dim(2)]);
    set(transv.h(12),'Xdata',[0 dim(1)],'Ydata',[1 1]*(dim(2)-P(2)+1));
   
    % colorbar
    %----------------------------------------------------------------------
    set(transv.h(14), 'Ydata',[mn mx], 'Cdata',(1:D)' + D);
    set(transv.h(13),'XTickLabel',[],'Ylim',[mn mx]);
    
else
    set(transv.h(2),'Cdata',rot90(spm_grid(T2)));
    set(get(transv.h(1),'Title'),'String',sprintf('z = %0.0fmm',xyzmm(3)));
    set(transv.h(3),'Xdata',[1 1]*P(1),'Ydata',[0 dim(2)]);
    set(transv.h(4),'Xdata',[0 dim(1)],'Ydata',[1 1]*(dim(2)-P(2)+1));

    % colorbar
    %----------------------------------------------------------------------
    set(transv.h(6), 'Ydata',[0 d], 'Cdata',(1:D)' + D);
    set(transv.h(5),'XTickLabel',[],'Ylim',[0 d]);
    
end


%-Reset pointer
%--------------------------------------------------------------------------
spm('Pointer','Arrow')
return;

%==========================================================================
% function clear_global(varargin)
%==========================================================================
function clear_global(varargin)
global transv
if isstruct(transv),
    for h = transv.h,
        if ishandle(h), set(h,'DeleteFcn',''); end;
    end
    for h = transv.h,
        if ishandle(h), delete(h); end;
    end
    transv = [];
    clear global transv;
end
return;
