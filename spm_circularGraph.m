function spm_circularGraph(A,varargin)
% Plot a circular graph to illustrate connections
% FORMAT spm_circularGraph(A,'PropertyName',propertyvalue,...)
% X      - symmetric (NxN) matrix of numeric or logical values
%
% Optional properties:
%   'Colormap'     - (Nx3) matrix of [r g b] triples
%   'Label'        - cell array of N strings
%
% A 'circular graph' is a visualization of a network of nodes and their
% connections. The nodes are laid out along a circle, and the connections
% are drawn within the circle.

% https://github.com/paul-kassebaum-mathworks/circularGraph

% Copyright (c) 2016, The MathWorks, Inc.
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
% * Redistributions of source code must retain the above copyright
%   notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
%   notice, this list of conditions and the following disclaimer in
%   the documentation and/or other materials provided with the distribution
% * Neither the name of the The MathWorks, Inc. nor the names
%   of its contributors may be used to endorse or promote products derived
%   from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


% Input parameters
%--------------------------------------------------------------------------
if ~(isnumeric(A) || islogical(A))
    error('Input A must be a symmetric matric of numeric or logical values.');
end

if strcmpi(spm_check_version,'matlab')
    ColorMap = parula(length(A));
else
    ColorMap = viridis(length(A));
end
Label = arrayfun(@(x) num2str(x), 1:length(A), 'UniformOutput',false);
for i=1:2:numel(varargin)
    switch lower(varargin{i})
        case 'colormap'
            ColorMap = varargin{i+1};
        case 'label'
            Label = varargin{i+1};
        otherwise
            error('Unknown option.');
    end
end
if length(ColorMap) ~= length(A) || length(Label) ~= length(A)
    error('Optional input length does not match input data.');
end
ax = gca;

% Draw the nodes
%--------------------------------------------------------------------------
t     = linspace(-pi,pi,length(A) + 1).'; % theta for each node
for i = 1:length(A)
    line(ax, cos(t(i)), sin(t(i)), 1,...
        'Marker', '.',...
        'MarkerSize', 64,...
        'Color', ColorMap(i,:));
    text(ax, cos(t(i)), sin(t(i)), 1,...
        Label{i},...
        'HorizontalAlignment', 'center');
end

% Find non-zero values of A and their indices
%--------------------------------------------------------------------------
A             = abs(A);
A             = A - diag(diag(A));
[row,col,v]   = find(A);

% Calculate line widths based on values of A (stored in v)
%--------------------------------------------------------------------------
minLineWidth  = 0;
lineWidthCoef = 4;
lineWidth     = v./max(v);
lineWidth     = lineWidthCoef*lineWidth + minLineWidth;

% Draw connections on the Poincare hyperbolic disk
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
            line(ax,...
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
            
            line(ax,...
                r*cos(theta)+x0,...
                r*sin(theta)+y0,...
                'LineWidth', lineWidth(i),...
                'Color', ColorMap(j,:));
        end
    end
end
