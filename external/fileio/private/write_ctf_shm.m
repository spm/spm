function [varargout] = write_ctf_shm(varargin)

% WRITE_CTF_SHM writes metainformation and data as a packet to shared memory.
% This function can be used for real-time processing of data while it is
% being acquired.
%
% Use as
%   write_ctf_shm(msgType, msgId, sampleNumber, numSamples, numChannels, data);

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: write_ctf_shm.m,v $
% Revision 1.1  2009/01/14 09:12:16  roboos
% The directory layout of fileio in cvs sofar did not include a
% private directory, but for the release of fileio all the low-level
% functions were moved to the private directory to make the distinction
% between the public API and the private low level functions. To fix
% this, I have created a private directory and moved all appropriate
% files from fileio to fileio/private.
%
% Revision 1.1  2007/08/01 09:38:28  roboos
% new implementation as mex file
%

error('could not locate MEX file');

