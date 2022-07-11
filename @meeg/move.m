function res = move(this, fname)
% Method for moving or changing name of data file
% FORMAT res = move(this, fname)
%
% fname can be
% - path\filename -> data moved and renamed
% - path          -> data moved only
% - filename      -> data renamed only
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


res = copy(this, fname);
delete(this);
