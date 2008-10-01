% BUFFER manages and accesses the realtime data acquisition buffer
% This function is implented as mex file.
%
% Use as
%   retval = buffer(cmd, detail, host, port)
%
% To read data from a buffer server over the network
%   hdr = buffer('get_hdr', [],     host, port)
%   dat = buffer('get_dat', datsel, host, port)
%   evt = buffer('get_evt', evtsel, host, port)
%
% The selection for data and events should be zero-offset and contain
%   datsel = [begsample endsample]
%   evtsel = [begevent  endevent]
%
% To write data to a buffer server over the network
%   buffer('put_hdr', hdr, host, port)
%   buffer('put_dat', dat, host, port)
%   buffer('put_evt', evt, host, port)
%
% To implement a local buffer server and have other clients
% connect to it at a specified network port
%   buffer('tcpserver', 'init', [], port)
%   buffer('tcpserver', 'exit', [], port)
%
% To implement a local acquisition client and have it buffer the
% data locally
%   buffer(acqclient, 'init')
%   buffer(acqclient, 'exit')
%
% To implement a local acquisition client and have it send the data
% through the network to a buffer server
%   buffer(acqclient, 'init', host, port)
%   buffer(acqclient, 'exit', host, port)
%
% Supported acquisition clients should be
%  'sinewave'
%  'ctf'
%  'brainamp'
%  'biosemi'
%  'tmsi'
%  'eldith'

% Copyright (C) 2008, Robert Oostenveld
%
% $Log: buffer.m,v $
% Revision 1.4  2008/07/08 20:31:34  roboos
% extended the list of acquisition threads in the mex file, also added all of them to thread stopping
%
% Revision 1.3  2008/05/22 09:23:09  roboos
% updated documentation
%
% Revision 1.2  2008/03/10 10:41:58  roboos
% updated documentation
%
% Revision 1.1  2008/03/09 22:53:04  roboos
% new function, only placeholder for the documentation since the actual implementation is in a mex file
%

error('could not locate MEX file');

