% FIELDTRIPDEFS is called at the begin of all FieldTrip functions and
% contains some defaults and path settings
%
% $Log: fieldtripdefs.m,v $
% Revision 1.6  2008/10/01 08:22:40  roboos
% added specest module
%
% Revision 1.5  2008/09/23 07:54:09  roboos
% added template directory, moved all layouts and cortical meshes there
%
% Revision 1.4  2008/09/23 07:40:10  roboos
% changed order of the path additioons, fixpath should be the first (i.e. lowest on the path)
%
% Revision 1.3  2008/09/23 07:31:13  roboos
% Moved all statfuns and trialfuns to their own directories, where they will be easier to find for the end-user. Also updated fieldtripdefs accordingly.
%
% Revision 1.2  2008/09/22 20:58:57  roboos
% added modules: fileio, forwinv, preproc
%
% Revision 1.1  2008/09/22 20:02:50  roboos
% initial version, contains fixpath
%

try
  % this is to help in converting the directory layout while at the same time maintainiong backward compatibility
  % NOTE THAT THIS ONE SHOULD BE ADDED THE FIRST TO ENSURE THAT IT IS THE LOWEST ON THE PATH
  hastoolbox('fixpath', 1, 1);
end

try
  % this contains layouts and cortical meshes
  hastoolbox('template', 1, 1);
end

try
  % this is used in statistics
  hastoolbox('statfun', 1, 1);
end

try
  % this is used in definetrial
  hastoolbox('trialfun', 1, 1);
end

try
  % numerous functions depend on this module
  hastoolbox('fileio', 1, 1);
end

try
  % numerous functions depend on this module
  hastoolbox('forwinv', 1, 1);
end

try
  % numerous functions depend on this module
  hastoolbox('preproc', 1, 1);
end

try
  % only very few functions depent on this module
  % it is not yet included in the FTP release version
  hastoolbox('specest', 1, 1);
end

