function [dat] = read_24bit(filename, offset, numwords)

% READ_24BIT read a stream of 24 bit values and converts them to doubles
% This function is designed for Biosemi BDF files and is implemented as mex
% file for efficiency.
%
% Use as
%   [dat] = read_24bit(filename, offset, numwords);
%
% See also READ_16BIT, READ_BIOSEMI_BDF

% Copyright (C) 2006, Robert Oostenveld
% F.C. Donders Ccentre for Cognitive Neuroimaging
% http://oase.uci.ru.nl/~roberto/
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: read_24bit.m,v $
% Revision 1.1  2006/02/01 08:34:10  roboos
% wrappers as placeholder for the help, the functionality is implemented as mex file
%

error('this function is implemented as mex file and is not available for this platform');
