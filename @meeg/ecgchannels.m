function ind = ecgchannels(this)
% Method for getting index vector of ecg channels
% FORMAT ind = ecgchannels(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Christophe Phillips & Stefan Kiebel
% $Id: ecgchannels.m 2668 2009-01-29 12:11:54Z christophe $

type = chantype(this);
ind = union(find(strcmpi('ECG', type)), find(strcmpi('EKG', type)));
ind = ind(:)'; % must be row to allow to use it as loop indices

