function [mri] = read_fcdc_mri(filename)

% READ_FCDC_MRI reads anatomical and functional MRI data from different
% file formats. The output data is structured in such a way that it is
% comparable to a FieldTrip source reconstruction.
%
% Use as
%   [mri] = read_fcdc_mri(filename)
%
% The output MRI may have a homogenous transformation matrix that converts
% the coordinates of each voxel (in xgrid/ygrid/zgrid) into head
% coordinates. 
%
% See also READ_FCDC_DATA, READ_FCDC_HEADER, READ_FCDC_EVENT

% Copyright (C) 2004-2009, Robert Oostenveld
%
% $Log: read_fcdc_mri.m,v $
% Revision 1.20  2009/02/17 11:33:01  roboos
% renamed read_fcdc_mri to read_mri and moved to fileio
% the old function remains as a wrapper
%

fieldtripdefs

% use the low-level reading function
[mri] = read_mri(filename);

