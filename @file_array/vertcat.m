function o = vertcat(varargin)
% Vertical concatenation of file_array objects.
% _______________________________________________________________________
% %W% John Ashburner %E%
o = cat(1,varargin{:});
return;
