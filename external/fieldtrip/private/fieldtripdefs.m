% FIELDTRIPDEFS is called at the begin of all FieldTrip functions and
% contains some defaults and path settings
%
% $Log: fieldtripdefs.m,v $
% Revision 1.1  2008/09/22 20:02:50  roboos
% initial version, contains fixpath
%

try
  % this is to help in converting the directory layout
  % while at the same time maintainiong backward compatibility
  hastoolbox('fixpath', 1, 1);
end

