function val = file2mat(a,varargin)
% Function for reading from file_array objects
% FORMAT val = file2mat(a,ind1,ind2,ind3,...)
% a      - file_array object
% indx   - indices for dimension x (int64)
% val    - the read values
%
% This function is normally called by file_array/subsref.
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: file2mat.m 7147 2017-08-03 14:07:01Z spm $

%-This is merely the help file for the compiled routine
error('file2mat.c not compiled - see Makefile');
