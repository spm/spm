function [out] = deblank2(in)

% DEBLANK2 removes all blanks from the beginning and end of a char array
%
% out = deblank2(in)

% Copyright (C) 2002-2007, Robert Oostenveld
%
% $Log: deblank2.m,v $
% Revision 1.3  2007/12/20 19:08:06  roboos
% changed some whitespace
%
% Revision 1.2  2003/03/17 10:37:28  roberto
% improved general help comments and added copyrights
%

out = fliplr(deblank(fliplr(deblank(in))));

