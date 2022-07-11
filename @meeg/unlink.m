function this = unlink(this)
% Unlinks the object from the data file 
% FORMAT this = unlink(this)   
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


this.data = [];
this      = check(this);
