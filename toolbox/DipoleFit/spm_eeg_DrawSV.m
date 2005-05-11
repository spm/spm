function [SV,h] = spm_eeg_DrawSV(SV,fig,name_el,label_t,film,paper)
%__________________________________________________________
%
% FORMAT [SV,h] = spm_eeg_DrawSV(SV,fig,name_el,label_t,film)
% 
% Function that displays the "Spline View" (SV) structure.
% The maps obtained are a simple flat disk representing the upper
% part of a sphere. Therefore, beware of the deformation when interpreting
% the maps !!!
% Maps can be drawn either as a single averaged-over-time map, as a series 
% of maps (one map for each time bin) or as a movie
% IN:
%   - SV      : "Spline Structure" to be displayed
%   - fig     : specify where the map is produced ([] in the spm graphics,
%               0 in a new figure window (default), or in figure 'fig')
%   - name_el : display the name of the electrodes (1) or not (0, default)
%   - label_t : Labels for the different plots
%   - film    : Create a movie from the time series (1) or not (0, default)
%   - paper   : Use larger font and lines (1) such as for paper figure,
%               or not (0, default)
%
% To re-display the movie, use 
%   n = 2 % Nr of repetitions (if <0, it goes back and forth)
%   fps = 10 % frames per second
%   movie(SV.hFigure,SV.M,n,fps,SV.pos)
%
% To generate an .avi movie file from the Matlab movie, use
%   movie2avi(SV.M,'filename.avi')
% For more options on quality/frames/etc use: help movie2avi
%
% All this has been adapted from routines written by Morgan Willis
% when he was working with Mick Rugg at the ICN, London, UK.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips,
% $Id: spm_eeg_DrawSV.m 143 2005-05-11 17:13:13Z christophe $

if nargin<6, paper = 0; end     % Normal font and line width
if nargin<5, film = 0; end      % No film display
if nargin<4, label_t = []; end  % No labels used
if nargin<3, name_el = 0; end   % No display electrodes name
if nargin<2, fig = []; end      % Use spm graph window (to create new window, 0)
if nargin<1, SV = spm_eeg_CreateSV; end % Create stuff to display

Nmaps = size(SV.interpolated,3);

% Check I have enough labels
if ~isempty(label_t)
    if length(label_t)==Nmaps
        if isnumeric(label_t)
            if size(label_t,1)==1
                label_t = num2str(label_t');
            else
                label_t = num2str(label_t);
            end
        elseif iscellstr(label_t)
            label_t = strvcat(label_t);
        else
            if ~ischar(label_t)
                warning('Could not figure out the label format, so they''re not displayed.');
                label_t = [];
            end
        end
    else
        % No label if not enough of them
        label_t = [];
    end
end

% WS = spm('WinScale');
if isempty(fig)
    tmp = spm_figure('FindWin');
    if isempty(tmp)
        hFig = spm_figure;
    else
        hFig = tmp;
    end
    spm_figure('ColorMap','jet')
elseif fig == 0
    hFig = figure;
    if film
        colormap('jet')
%         colormap('gray')
    end
else
    hFig = figure(fig);
end
figure(hFig);

if film
    bord_sp = .1;
    betw_sp = .03;
    Nc = 1; Nr = 1;
	hori_sp = 1-2*bord_sp;
	vert_sp = 1-2*bord_sp;
else
    bord_sp = .04;
    betw_sp = .03;
	maxNrColumn = ceil(sqrt(Nmaps));
	Nc = min(maxNrColumn,Nmaps);
	Nr = ceil(Nmaps/Nc);
	hori_sp = (1-2*bord_sp-(Nc-1)*betw_sp)/Nc;
	vert_sp = (1-2*bord_sp-(Nr-1)*betw_sp)/Nr;
end

if paper
    size2use = [2,16,8]; % Line width, font size and marker size
else
    size2use = [1,10,4]; % Line width, font size and marker size
end


% Prepare eletrodes location on map
X = SV.spline.X;
[xf,yf] = CartToFlat( X(1,:), X(2,:), X(3,:) );

% Loop over the images
for ii=1:Nmaps
    if ~film
        r_ii = floor((ii-1)/Nc)+1;
        c_ii = mod(ii-1,Nc)+1;
    elseif film & ii==1
        r_ii = 1 ;
        c_ii = 1 ;
    end
    
    if (film & ii==1) | ~film
        temp_axis = axes('Position',[bord_sp+(c_ii-1)*(hori_sp+betw_sp) ...
                                    1-bord_sp-vert_sp-(r_ii-1)*(vert_sp+betw_sp) hori_sp vert_sp]);
    end
    if film & ii==1
        hAx = temp_axis;
    elseif ~film
        hAx(ii) = temp_axis;
    end
    if (film & ii==1) | ~film
        % Draw interpolated values
        set( temp_axis, 'XLim', [-SV.maxphi SV.maxphi], 'YLim', [-SV.maxphi SV.maxphi] );
        axis image
        s = size( SV.interpolated(:,:,ii) );
        xg = -SV.maxphi  :  2*SV.maxphi / (s(2)-1)  :  SV.maxphi;
        yg = -SV.maxphi  :  2*SV.maxphi / (s(1)-1)  :  SV.maxphi;
        hSurf = surface( xg, -yg, zeros(s), SV.interpolated(:,:,ii), 'FaceColor', 'interp', ...
            'EdgeColor', 'none', 'HitTest', 'off'  );
        axis off
        
        % Draw circles
        rectangle( 'Position', [ -pi/4 -pi/4 pi/2 pi/2 ],         ...
            'Curvature', [1 1], 'LineWidth', size2use(1), 'HitTest', 'off'  );
        rectangle( 'Position', [ -pi/2 -pi/2 pi pi ],             ...
            'Curvature', [1 1], 'LineWidth', size2use(1), 'HitTest', 'off' );
        rectangle( 'Position', [ -SV.maxphi -SV.maxphi SV.maxphi*2 SV.maxphi*2 ], ...
            'Curvature', [1 1], 'LineWidth', size2use(1), 'HitTest', 'off' );
        if name_el
            text(0,-pi/4-.05,'45^{o}','FontSize',size2use(2),'Interpreter','tex');
            text(0,-pi/2-.05,'90^{o}','FontSize',size2use(2),'Interpreter','tex');
            text(0,-SV.maxphi-.06,[num2str(SV.maxphi/pi*180),'^{o}'], ...
                    'FontSize',size2use(2),'Interpreter','tex');
        end
        
        % Draw lines
        line( [ 0 0 ], [-SV.maxphi SV.maxphi], 'Color', [ 0 0 0 ], ...
            'LineWidth', size2use(1), 'HitTest', 'off' );
        line( [ -SV.maxphi SV.maxphi ], [0 0], 'Color', [ 0 0 0 ], ...
            'LineWidth', size2use(1), 'HitTest', 'off');
        
        % Draw electrodes
        hold on;
        plot( xf, yf, 'ow', 'LineWidth', 1, 'MarkerFaceColor', 'k', ...
            'MarkerSize', size2use(3), 'HitTest', 'off' );
        hold off;
        
        % Display electrodes' name
        if name_el
            text(xf,yf,SV.el_name);
        end
    
        %Colorbar
        caxis(SV.range);
        if film
            hCol = colorbar;
        else
            hCol(ii) = colorbar;
        end
        if paper
            set(hCol(end),'Fontsize',size2use(2),'Linewidth',2)
        end
	
        % Add label
        if ~isempty(label_t)
            hTitl = title(label_t(ii,:));
        end
	
    elseif film & ii>1
        set(hSurf,'CData',SV.interpolated(:,:,ii));
        if ~isempty(label_t)
            set(hTitl,'String',label_t(ii,:));
        end
    end    
    % If a film is created
    if film
        drawnow;
        SV.M(:,ii) = getframe(hFig);
    end
end

if film
    delete(hAx);
    movie(hFig,SV.M,film)
end

SV.hFigure = hFig;
SV.hAxis = hAx;
SV.hColourBar = hCol;
h = struct('hFig',hFig,'hAx',hAx,'hCol',hCol);
SV.pos = [bord_sp bord_sp 1 1];

%_________________________________________________________________________________
%
% SUB_FUNCTIONS :
%_________________________________________________________________________________
function [ xf, yf ] = CartToFlat( x, y, z )
% Convert Cartesian coordinates on surface of unit sphere to 2D Cartesian 
% Flat Map coordinates


rh = sqrt( x.*x + y.*y );

zero = find( rh == 0 );
rh( zero ) = 1e-20;

costheta = x ./ rh;
sintheta = y ./ rh;
    
%phi   = atan(  sqrt( x.*x + y.*y )  ./  z  );

% Modification by chrisp@fil, 02/02/24,
% to take into account electrodes with phi<-90 or phi>90
%---------------------------------------------------------------------
lz = find(z<0) ;
% phi   = atan(  sqrt( x.*x + y.*y )  ./  z  );
phi   = atan2( sqrt( x.*x + y.*y ) , z );
% phi(lz) = pi + atan(  sqrt( x(lz).*x(lz) + y(lz).*y(lz) )  ./  z(lz)  );

xf = costheta .* phi;
yf = sintheta .* phi;



