function plot_sens(sens, varargin)

% PLOT_SENS plots the position of the channels in the EEG or MEG sensor array
%
% Use as
%   plot_sens(sens, ...)
% where optional input arguments should come in key-value pairs and can
% include
%   'style'    plotting style for the points representing the channels, see plot3 (default = 'k.')

% Copyright (C) 2009, Robert Oostenveld
%
% $Log: plot_sens.m,v $
% Revision 1.3  2009/04/14 19:48:28  roboos
% added keyvalcheck
%
% Revision 1.2  2009/04/08 06:35:05  roboos
% first implementation, covers the use in headmodelplot
%

% get the optional input arguments
keyvalcheck(varargin, 'optional', {'style'});
style = keyval('style', varargin); if isempty(style), style = 'k.'; end

% determine the position of each channel, which is for example the mean of
% two bipolar electrodes, or the bottom coil of a axial gradiometer
[chan.pnt, chan.label] = channelposition(sens);

plot3(chan.pnt(:,1), chan.pnt(:,2), chan.pnt(:,3), style);

