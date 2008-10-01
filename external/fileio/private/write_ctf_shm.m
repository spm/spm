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
% Revision 1.1  2007/08/01 09:38:28  roboos
% new implementation as mex file
%

error('could not locate MEX file');

