function res = fiducials(this, newfiducials)
% Method for getting/setting the fiducials field
% FORMAT res = fiducials(this, fiducials)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id $

switch nargin
    case 1
         res = this.fiducials;
    case 2
         this.fiducials = newfiducials;
         res = this;
end
