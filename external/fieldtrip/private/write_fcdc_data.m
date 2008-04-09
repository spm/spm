function write_fcdc_data(varargin)

% WRITE_FCDC_DATA exports electrophysiological data to a file
%
% Use as
%   write_fcdc_data(filename, dat, ...)
%
% The specified filename can already contain the filename extention,
% but that is not required since it will be added automatically.
%
% Additional options should be specified in key-value pairs and can be
%   'dataformat'     string, see below
%   'header'         header structure, see READ_FCDC_HEADER
%   'chanindx'       1xN array, selection of channels from the header
%
% The supported dataformats are
%   brainvision_eeg
%   neuralynx_ncs
%   plexon_nex
%   fcdc_matbin
%   matlab

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: write_fcdc_data.m,v $
% Revision 1.5  2007/06/13 08:14:31  roboos
% updated documentation
%
% Revision 1.4  2007/06/13 06:39:37  roboos
% moved the content of the write_fcdc_data to the low-level write_data function
% updated the help
%
% Revision 1.3  2007/03/20 17:05:40  roboos
% removed unused subfunction for fixing the file extension
%
% Revision 1.2  2007/03/20 17:04:27  roboos
% reimplemented plexon_nex with a nex-structure
% ensure that the extension is correct
% added plain matlab as output format
% moved the scaling and conversion to int16 into the low-level function (for neuralynx)
% updated the documentation
%
% Revision 1.1  2007/02/21 09:50:56  roboos
% initial implementation with support for fcdc_matbin, neuralynx_ncs, plexon_nex
%

% use the low-level writing function
write_data(varargin{:});
