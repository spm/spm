function [varargout] = read_ctf_shm(varargin);

% READ_CTF_SHM reads metainformation or selected blocks of data from
% shared memory. This function can be used for real-time processing of
% data while it is being acquired.
%
% Use as
%   [msgType msgId sampleNumber numSamples numChannels] = read_ctf_shm;
% or
%   [data] = read_ctf_shm(msgNumber);
%   [data] = read_ctf_shm(msgNumber, numValues);

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: read_ctf_shm.m,v $
% Revision 1.1  2009/01/14 09:12:15  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.2  2007/07/30 11:51:55  roboos
% updated documentation
%
% Revision 1.1  2007/07/24 11:17:03  roboos
% first implementation, tested and works on odin
%

error('could not locate MEX file');

