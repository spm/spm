function mat2file(a,val,varargin)
% Function for writing to file_array objects.
% FORMAT mat2file(a,val,ind1,ind2,ind3,...)
% a      - file_array object
% val    - values to write
% indx   - indices for dimension x (int32)
%
% This function is normally called by file_array/subsasgn
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: mat2file.m 201 2005-07-20 10:57:13Z guillaume $

%-This is merely the help file for the compiled routine
error('mat2file.c not compiled - see Makefile');
