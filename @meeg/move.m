function res = move(this, fname)
% Method for moving or changing name of data file
% FORMAT res = move(this, fname)
%
% fname can be
% - path\filename -> data moved and renamed
% - path          -> data moved only
% - filename      -> data renamed only
%__________________________________________________________________________
% Copyright (C) 2011-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: move.m 5025 2012-10-31 14:44:13Z vladimir $

res = copy(this, fname);
delete(this);
