function [p,f] = spm_eeg_inv_displScEl(varargin)

% p = spm_eeg_inv_displScEl(scalp,electr) or spm_eeg_inv_displScEl(model)
%
% Display a tessalated surface, in color, with the electrodes on top.
% 'spm_eeg_inv_displScEl' returns handle to patch of the scalp
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips,
% $Id$

if nargin < 2
    if nargin == 0
        Pmod = spm_select(1,'^model.*e\d.*\.mat','Select model to display');
        load(Pmod)
    else
        model = varargin{1};
    end
    if isfield(model,'head') & isfield(model,'electrodes')
        scalp  = model.head(length(model.head));
        electr = model.electrodes;
    else
        error('Wrong input arguments')
    end
elseif nargin == 2
    scalp  = varargin{1};
    electr = varargin{2};
else
    error('Wrong input arguments')
end

Nvert = scalp.nr(1) ;
Ntri  = scalp.nr(2) ;
vert  = scalp.XYZmm';
tri   = scalp.tri'  ;

% Generate the electrodes.
[x_sph,y_sph,z_sph] = sphere ;
x_sph = x_sph*3 ; y_sph = y_sph*3 ; z_sph = z_sph*3 ; 
c_sph = ones(size(x_sph))*5 ;

colo_skin = [1 .7 .55] ;
c    = ones(Ntri,1)*colo_skin ;
hFig = spm_figure('FindWin','Graphics');
if isempty(hFig)
    hFig = spm_figure;
else
    spm_figure('Clear',hFig);
    spm_figure('ColorMap','gray');
end

figure(hFig);
% set(hFig,'Render','OpenGL')
set(hFig,'Render','zbuffer')
p = patch('Vertices',vert,'Faces',tri,'FaceColor',colo_skin) ;
set(p,'FaceLighting','phong','SpecularStrength',.1,...
	'AmbientStrength',.45,'EdgeColor','none'); %,'FaceColor','interp') ;

el_coord = electr.XYZmm ;

hold on
for ii=1:electr.nr
	s = surf(x_sph+el_coord(1,ii),y_sph+el_coord(2,ii),...
		z_sph+el_coord(3,ii),c_sph) ;
	set(s,'EdgeColor','none') ;
end
caxis([0 5])

view(3)
axis equal
axis vis3d
axis off

view(135,15)

h1 = light('Position',[1 3 3]) ;
h2 = light('Position',[1 -1 -1]) ;
h3 = light('Position',[-3 -1 0]) ;

rotate3d on
