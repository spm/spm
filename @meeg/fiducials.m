function res = fiducials(this, newfiducials)
% Method for getting/setting the fiducials field
% FORMAT res = fiducials(this, fiducials)
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: fiducials.m 5025 2012-10-31 14:44:13Z vladimir $

switch nargin
    case 1
         res = this.fiducials;
    case 2
         this.fiducials = newfiducials;
         res = this;
end
