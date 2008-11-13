function [varargout] = deepcopy(varargin)

% DEEPCOPY makes a deep copy of an array, and returns a pointer to
% the copy. A deep copy refers to a copy in which all levels of data
% are copied. For example, a deep copy of a cell array copies each
% cell, and the contents of the each cell (if any), and so on.
%
% Example
%   clear a b c
%   a = 1;
%   b = a;            % this is a regular copy
%   c = deepcopy(a);  % this is a deep copy
%   increment(a);     % increment the value of a with one, using pass by reference
%   disp(a);
%   disp(b);
%   disp(c);
%  

error('Could not locate mex file.');

