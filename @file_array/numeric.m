function out = numeric(fa)
% Convert to numeric form
% FORMAT numeric(fa)
% fa - a file_array
% _______________________________________________________________________
% %W% John Ashburner %E%

[vo{1:ndims(fa)}] = deal(':');
out = subsref(fa,struct('type','()','subs',{vo}));

