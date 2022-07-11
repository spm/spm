function res = fiducials(this, newfiducials)
% Method for getting/setting the fiducials field
% FORMAT res = fiducials(this, fiducials)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


switch nargin
    case 1
         res = this.fiducials;
    case 2
         this.fiducials = ft_struct2double(fixpnt(newfiducials));
         res = this;
end
