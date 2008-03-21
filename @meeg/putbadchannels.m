function obj = putbadchannels(obj, ind, flag)
% Method for chaing a rejection flag
% FORMAT res = putbadchannels(obj, ind, flag)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

[obj.channels(ind).bad] = deal(flag);