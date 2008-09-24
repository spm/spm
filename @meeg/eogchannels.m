function ind = eogchannels(this)
% Method for getting index vector of eog channels
% FORMAT ind = eogchannels(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: eogchannels.m 2164 2008-09-24 11:48:57Z stefan $

type = chantype(this);
ind = union(find(strcmpi('EOG', type)), find(strcmpi('HEOG', type)), find(strcmpi('VEOG', type)));
ind = ind(:)'; % must be row to allow to use it as loop indices

