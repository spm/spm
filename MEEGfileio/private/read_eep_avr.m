function [dat] = read_eep_avr(fn);

% READ_EEP_AVR reads averaged EEG data from an EEProbe *.avr file
% and returns a structure containing the header and data information.
%
% eeg = read_eep_avr(filename)
%
% See also READ_EEP_CNT, READ_EEP_TRG, READ_EEP_REJ

% Copyright (C) 2002, Robert Oostenveld
%
% $Log: read_eep_avr.m,v $
% Revision 1.3  2003/03/13 14:25:45  roberto
% converted from DOS to UNIX
%
% Revision 1.2  2003/03/11 15:24:51  roberto
% updated help and copyrights
%

error('could not locate mex file');
