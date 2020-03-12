function spm_circularGraph(A,varargin)
% spm_circularGraph: Plots a circular graph to illustrate connections
%
%% Syntax
% spm_circularGraph(X)
% spm_circularGraph(X,'PropertyName',propertyvalue,...)
% h = spm_circularGraph(...)
%
%% Description
% A 'circular graph' is a visualization of a network of nodes and their
% connections. The nodes are laid out along a circle, and the connections
% are drawn within the circle.
% Required input arguments.
% X  (N,N) : A symmetric matrix of numeric or logical values.
%
% Optional properties.
% Colormap : A N by 3 matrix of [r g b] triples, where N is the
%            length(adjacenyMatrix).
%%
% Copyright 2016 The MathWorks, Inc.

% Constructor
%--------------------------------------------------------------------------
p = inputParser;

defaultColorMap = parula(length(A));
addRequired(p,'A',@(x)(isnumeric(x) || islogical(x)));
addParameter(p,'ColorMap',defaultColorMap,@(colormap)length(colormap) == length(A));

parse(p,A,varargin{:});
ColorMap = p.Results.ColorMap;

% Draw the nodes
%--------------------------------------------------------------------------
t     = linspace(-pi,pi,length(A) + 1).'; % theta for each node
for i = 1:length(A)
    line(cos(t(i)),sin(t(i)),1,...
        'Marker','.','MarkerSize',64,...
        'Color',ColorMap(i,:))
        text(cos(t(i)),sin(t(i)),1,num2str(i),'HorizontalAlignment','center')
end

% Find non-zero values of s and their indices
%--------------------------------------------------------------------------
A             = abs(A);
A             = A - diag(diag(A));
[row,col,v]   = find(A);

% Calculate line widths based on values of s (stored in v)
%--------------------------------------------------------------------------
minLineWidth  = 0;
lineWidthCoef = 4;
lineWidth     = v./max(v);
lineWidth     = lineWidthCoef*lineWidth + minLineWidth;

% Draw connections on the Poincare hyperbolic disk.
%==========================================================================
% Equation of the circles on the disk:
% x^2 + y^2
% + 2*(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1))*x
% - 2*(u(1)-v(1))/(u(1)*v(2)-u(2)*v(1))*y + 1 = 0,
% where u and v are points on the boundary.
%
% Standard form of equation of a circle
% (x - x0)^2 + (y - y0)^2 = r^2
%
% Therefore we can identify
% x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
% y0 = (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
% r^2 = x0^2 + y0^2 - 1

for i = 1:length(v)
    if row(i) ~= col(i) && lineWidth(i) > 1
        
        % color index (of the stongest efferent connection)
        %------------------------------------------------------------------
        if A(row(i),col(i)) > A(col(i),row(i))
            j = col(i);
        else
            j = row(i);
        end
        
        
        % draw lines
        %------------------------------------------------------------------
        if abs(row(i) - col(i)) - length(A)/2 == 0
            
            % points are diametric, so draw a straight line
            %--------------------------------------------------------------
            u = [cos(t(row(i)));sin(t(row(i)))];
            v = [cos(t(col(i)));sin(t(col(i)))];
            line(...
                [u(1);v(1)],...
                [u(2);v(2)],...
                'LineWidth', lineWidth(i),...
                'Color', ColorMap(j,:));
        else
            
            % points are not diametric, so draw an arc
            %--------------------------------------------------------------
            u  = [cos(t(row(i)));sin(t(row(i)))];
            v  = [cos(t(col(i)));sin(t(col(i)))];
            x0 = -(u(2)-v(2))/(u(1)*v(2)-u(2)*v(1));
            y0 =  (u(1)-v(1))/(u(1)*v(2)-u(2)*v(1));
            r  = sqrt(x0^2 + y0^2 - 1);
            thetaLim(1) = atan2(u(2)-y0,u(1)-x0);
            thetaLim(2) = atan2(v(2)-y0,v(1)-x0);
            
            % ensure the arc is within the unit disk
            %--------------------------------------------------------------
            if u(1) >= 0 && v(1) >= 0
                
                theta = [linspace(max(thetaLim),pi,50),...
                    linspace(-pi,min(thetaLim),50)].';
            else
                theta = linspace(thetaLim(1),thetaLim(2)).';
            end
            
            line(...
                r*cos(theta)+x0,...
                r*sin(theta)+y0,...
                'LineWidth', lineWidth(i),...
                'Color', ColorMap(j,:));
        end
    end
end
