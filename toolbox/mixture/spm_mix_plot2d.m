function [] = spm_mix_plot2d (mix,area,nContLines,LineType,min_p,max_p)
% Plot component density contours for 2D mixture model
% FORMAT [] = spm_mix_plot2d (mix,area,nContLines,LineType,min_p,max_p)
%
% mix           Mixture model data structure
% area          [xmin,xmax,ymin,ymax]
% nContLines    Number of contour lines; default=10
% LineType      Plot line type; default='-'
% min_p/max_p   Values of min and max probability contours
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny 
% $Id: spm_mix_plot2d.m 5647 2013-09-20 13:03:44Z ged $

if nargin < 2 | isempty(area), area=[0 10 0 10]; end
if nargin < 3 | isempty(nContLines), nContLines=10; end
if nargin < 4 | isempty(LineType), LineType='-'; end
if nargin < 5 | isempty(min_p), min_p=0.005; end
if nargin < 6 | isempty(max_p), max_p=0.1; end

xmin=area(1);
xmax=area(2);
ymin=area(3);
ymax=area(4);

% Number of data points per dimension
d=100; 

% Generate X   
dx1=(xmax-xmin)/d;
dx2=(ymax-ymin)/d;

x1=[xmin:dx1:xmax];
y1=[ymin:dx2:ymax];
[g1,g2]=meshgrid(x1,y1);

xplot = [reshape(g1,(d+1)^2,1), reshape(g2,(d+1)^2,1)];

held = ishold; cla; hold on

% Plot contours for each component
for j=1:mix.m,
    y = spm_MNpdf(mix.state(j).m, mix.state(j).C, xplot);
    
    % Plot proby contours
    yplot = reshape(y,d+1,d+1);
    dp=(max_p-min_p)/nContLines;
    clevels=[min_p:dp:max_p]*max(y);        
    
    clevels=clevels(1:nContLines);
    if nContLines==1
        contour(g1,g2,yplot,[clevels, clevels],LineType);
    else
        contour(g1,g2,yplot,clevels,LineType);
    end
end

if ~held, hold off, end
