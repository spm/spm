function this = rmfield(this, fields)
% Method for removing an object field
% FORMAT this = rmfield(this, fields)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


 this.other = rmfield(this.other, fields);
