function this = blank(this, fnamedat)
% Creates a blank datafile matching in the header in dimensions 
% Will not erase existing datafile it it's there
% FORMAT this = blank(this) 
%   Will create the datafile using fname and path     
% FORMAT this = blank(this, fnamedat) 
%   Will create the datafile using the provided name and path
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


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

initialise(this.data);

this = check(this);
