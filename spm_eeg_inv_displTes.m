function [p,hFig] = spm_eeg_inv_displTes(tsurf,c)
% FORMAT [p,f] = spm_eeg_inv_displTes(tsurf,c)
%
% Display a tessalated surface, with color c if provided : 
% if c is not provided, it uses 1 color for all surface.
% 'spm_eeg_inv_displTes' returns handle to patch & figure
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips,
% $Id$

% mono_c = 1/0 uses monocolor or provided c
% tri_c = 1/0 , color specified on triangles or vertices

if nargin==0
    Pmod = spm_select(1,'^model.*\.mat$','Select model to display');
    load(Pmod)
    if length(model.head)>1
        % Select surface to display
        list = '';
        for ii=1:length(model.head)
            list = strvcat(list,model.head(ii).info.str);
        end
        tsurf_i = spm_input('Display which surface','+1','m',list);
        tsurf = model.head(tsurf_i);
    else
        tsurf = model.head(1);
    end
end


Nvert = tsurf.nr(1) ;
Ntri = tsurf.nr(2) ;
try
    vert = tsurf.M*[tsurf.XYZvx ; ones(1,Nvert)]; vert = vert(1:3,:)';
catch
    vert = tsurf.XYZmm';
end
tri = tsurf.tri' ;

%colo_skin = [1 .7 .55] ;
%colo_skin = [1 .5 .45] ;
colo_skin = [1 0 0] ;

mono_c = 0 ;
if nargin<2
	c = ones(Ntri,1)*colo_skin ;
	mono_c = 1 ;
	tri_c = 1 ;
elseif nargin==2
	if (size(c,1)==1) & ((size(c,2)==1)|(size(c,2)==3))
		mono_c = 1 ;
	elseif size(c,1)==tsurf.nr(2)
		tri_c = 1 ;
		c3 = apply_colormap(c) ;
		c = c3 ;
	elseif size(c,1)==tsurf.nr(1)
		tri_c = 0 ;
		c3 = apply_colormap(c) ;
		c = c3 ;
	else
		error('Wrong colour specification') ;
	end
end

hFig = spm_figure('FindWin');
if isempty(hFig)
    hFig = figure;
else
    spm_figure('Clear',hFig);
    spm_figure('ColorMap','jet');
end

figure(hFig) ;
% set(f,'Render','OpenGL')
set(hFig,'Render','zbuffer')
if mono_c
	p = patch('Vertices',vert,'Faces',tri,'FaceColor',colo_skin) ;
%	p = patch('Vertices',vert,'Faces',tri,'FaceColor',[ 0.3 0.8 0.8 ]) ;
%	p = patch('Vertices',vert,'Faces',tri,'FaceVertexCData',c) ;
%	set(p,'FaceLighting','phong','SpecularStrength',.1,...
%		'AmbientStrength',.45,'EdgeColor','none') ;
		% ,'FaceColor','interp') ;
%	set(p,'FaceLighting','phong','SpecularStrength',.001,...
%		'AmbientStrength',.2,'EdgeColor','none',...
%		'DiffuseStrength',.6,'SpecularExponent',20) ;
		% ,'FaceColor','interp') ;
else
	if tri_c
		p = patch('Vertices',vert,'Faces',tri,'FaceVertexCData',c)
		set(p,'FaceLighting','phong','SpecularStrength',.1,...
			'AmbientStrength',.45,'EdgeColor','none',...
			'FaceColor','interp') ;
	else
		p = patch('Vertices',vert,'Faces',tri,'FaceVertexCData',c)
		set(p,'FaceLighting','phong','SpecularStrength',.1,...
			'AmbientStrength',.45,'EdgeColor','none',...
			'FaceColor','interp') ;
	end
end

view(3)
axis equal
%axis([1 181 1 217 1 181])
axis vis3d
rotate3d on
axis off

view(135,15)

%h1 = light('Position',[1 3 3]) ;
%h2 = light('Position',[1 -1 -1]) ;
%h3 = light('Position',[-3 -1 0]) ;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function c3 = apply_colormap(c)

% Scale between the min and max of the values in c
mc = min(c) ; Mc = max(c) ; dc = Mc-mc ;

% Scale symetrically around 0 the values of c
%mc = -max(abs(c)) ; dc = -2*mc ;

% Colormap : white-blue-black-red-white
%colormap(hot(128)), tmp_scale = colormap ; tmp2 = fliplr(flipud(tmp_scale)) ; col_scale = [tmp2 ; tmp_scale] ;

%colormap(jet(256)), col_scale = colormap ;
%colormap(gray(256)), col_scale = colormap ;
colormap(cool(256)), col_scale = colormap ;


colormap(col_scale)

ind_c = round((c-mc)/dc*255+1) ;

c3 = col_scale(ind_c,:) ;
