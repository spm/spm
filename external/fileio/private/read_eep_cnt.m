function [eeg] = read_eep_cnt(varargin);

% READ_EEP_CNT reads continuous EEG data from an EEProbe *.cnt file
% and returns a structure containing the header and data information.
%
% eeg = read_eep_cnt(filename, sample1, sample2)
%
% where sampel1 and sample2 are the begin and end sample of the data
% to be read.
%
% See also READ_EEP_TRG, READ_EEP_REJ, READ_EEP_AVR

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: read_eep_cnt.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.4  2004/03/29 15:18:51  roberto
% minor changes to the naming of input arguments
%
% Revision 1.3  2003/03/13 14:25:45  roberto
% converted from DOS to UNIX
%
% Revision 1.2  2003/03/11 15:24:51  roberto
% updated help and copyrights
%

error('could not locate mex file');
