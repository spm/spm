function out = isempty(this)
% True if the object is empty 
% FORMAT out = isempty(this)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


out = all(size(this)==0);
