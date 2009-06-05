function [x, y] = select_point2d(varargin)

% SELECT_POINT2d helper function for selecting a rectangular region
% in the current figure using the mouse.
%
% Use as
%   [x, y] = select_point2d(...)
%
% It returns a list of  [x  y] coordinates
% of the points in the selected region.
%
% Optional input arguments should come in key-value pairs and can include
%   'multiple'   true/false, make multiple selections by dragging, clicking
%                in one will finalize the selection (default = false)

% $Log: select_point2d.m,v $
% Revision 1.3  2009/06/04 10:51:09  roboos
% only whitespace
%
% Revision 1.2  2009/06/03 11:23:46  crimic
% first implementation
%

[box_x box_y] = select_box;

x = sort([box_x(1) box_x(2)]);
y = sort([box_y(1) box_y(2)]);  
