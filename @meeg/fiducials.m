function res = fiducials(this, newfiducials)
% Method for getting/setting the fiducials field
% FORMAT res = fiducials(this, fiducials)
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: fiducials.m 6622 2015-12-03 11:54:13Z vladimir $

switch nargin
    case 1
         res = this.fiducials;
    case 2
         this.fiducials = ft_struct2double(fixpnt(newfiducials));
         res = this;
end
