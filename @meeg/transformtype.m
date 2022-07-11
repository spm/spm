function res = transformtype(this, newtype)
% Method for getting/setting type of transform
% FORMAT res = transformtype(this, name)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if nargin == 1
    res = this.transform.ID;
else
    if strncmpi(this, 'TF', 2) && length(size(this))~=4
        error('TF transformtype can only be assigned to 4D dataset');
    end
    
    this.transform.ID = newtype;
    res = this;
end
