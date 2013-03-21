function this = blank(this, fnamedat)
% Creates a blank datafile matching in the header in dimensions 
% Will not erase existing datafile it it's there
% FORMAT this = blank(this) 
%   Will create the datafile using fname and path     
% FORMAT this = blank(this, fnamedat) 
%   Will create the datafile using the provided name and path
% _________________________________________________________________________
% Copyright (C) 2011-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: blank.m 5342 2013-03-21 15:50:05Z vladimir $

if nargin == 1
    [p, f] = fileparts(fullfile(this));
    fnamedat = fullfile(p, [f '.dat']);
else
   [p, f, x] = fileparts(fnamedat);
   if isempty(p)
       p = path(this);
   end
   if isempty(x)
       x = '.dat';
   end
   fnamedat = fullfile(p, [f x]); 
end

if exist(fnamedat, 'file')
    error('Data file exists. Use rmdata to delete.')
end

if isempty(this)
   error('All header dimensions should be >0');
end

this.data = file_array(fnamedat, size(this), 'float32-le');

endind = num2cell(size(this));
this.data(endind{:}) = 0;

this = check(this);
