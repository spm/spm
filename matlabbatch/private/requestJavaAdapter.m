function result = requestJavaAdapter(object)
%REQUESTJAVAADAPTER Support function for GUIDE

%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 1915 $
 
%%%%%%%  CAUTION                                           %%%%%%%%
%%%%%%%  This file is duplicated in both uitools and guide %%%%%%%%%
%%%%%%%  %TODO - determine if this functionality can be    %%%%%%%%%
%%%%%%%  broken out or replaced so this file doesn't exist %%%%%%%%%
%%%%%%%  in two places                                     %%%%%%%%%
 
  len = length(object);
  if len == 1
    if ishandle(object) && ~isjava(object)
      result = java(handle(object));
    else
      error('MATLAB:requestJavaAdapter:InvalidInput', '''requestJavaAdapter'' argument must be a handle list.');
    end
  else
    if ~isempty(object) && all(ishandle(object) & ~isjava(object))
      result = cell(len, 1);
      for i = 1:len
        result{i} = java(handle(object(i)));
      end
    else
      error('MATLAB:requestJavaAdapter:InvalidArgument', '''requestJavaAdapter'' argument must be a handle list.');
    end
  end
