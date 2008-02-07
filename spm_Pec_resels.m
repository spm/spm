function R = spm_Pec_resels(L,W)
% Returns the resel count for a point-list of voxels
% FORMAT R = spm_Pec_resels(L,W)
% L   - point list of voxels {in voxels}
% W   - smoothness of the component fields {FWHM in voxels}
% R   - vector of RESEL counts
%___________________________________________________________________________
%
% Reference : Worsley KJ et al 1996, Hum Brain Mapp. 4:58-73
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_Pec_resels.m 1143 2008-02-07 19:33:33Z spm $


% Resel Counts
%---------------------------------------------------------------------------
L   = round(L);
Ex  = 0;
Ey  = 0;
Ez  = 0;
Fxy = 0;
Fxz = 0;
Fyz = 0;
C   = 0;
x   = L(1,:);
y   = L(2,:);
z   = L(3,:);
P   = length(x);
R   = zeros(1,4);

% characterize voxel space
%---------------------------------------------------------------------------
for i = 1:P
  d = any(~any([x - L(1,i) - 1;y - L(2,i) - 0;z - L(3,i) - 0]));
  if d
    Ex = Ex + 1;
    d  =     any(~any([x - L(1,i) - 0;y - L(2,i) - 1;z - L(3,i) - 0]));
    d  = d & any(~any([x - L(1,i) - 1;y - L(2,i) - 1;z - L(3,i) - 0]));
    if d
    Fxy = Fxy + 1;
    d   =     any(~any([x - L(1,i) - 0;y - L(2,i) - 0;z - L(3,i) - 1]));
    d   = d & any(~any([x - L(1,i) - 1;y - L(2,i) - 0;z - L(3,i) - 1]));
    d   = d & any(~any([x - L(1,i) - 1;y - L(2,i) - 1;z - L(3,i) - 1]));
    d   = d & any(~any([x - L(1,i) - 0;y - L(2,i) - 1;z - L(3,i) - 1]));
    if d
        C = C + 1;
    end
    end
    d  =     any(~any([x - L(1,i) - 0;y - L(2,i) - 0;z - L(3,i) - 1]));
    d  = d & any(~any([x - L(1,i) - 1;y - L(2,i) - 0;z - L(3,i) - 1]));
    if d
    Fxz = Fxz + 1;
    end
  end
  d = any(~any([x - L(1,i) - 0;y - L(2,i) - 1;z - L(3,i) - 0]));
  if d
    Ey = Ey + 1;
    d  =     any(~any([x - L(1,i) - 0;y - L(2,i) - 0;z - L(3,i) - 1]));
    d  = d & any(~any([x - L(1,i) - 0;y - L(2,i) - 1;z - L(3,i) - 1]));
    if d
    Fyz = Fyz + 1;
    end
  end
  d = any(~any([x - L(1,i) - 0;y - L(2,i) - 0;z - L(3,i) - 1]));
  if d
    Ez = Ez + 1;
  end
  x(1) = [];
  y(1) = [];
  z(1) = [];
end

% Resel counts
%---------------------------------------------------------------------------
r     = 1./W(:);

R(1)  = P - (Ex + Ey + Ez) + (Fyz + Fxz + Fxy) - C;
R(2)  = (Ex - Fxy - Fxz + C)*r(1) + (Ey - Fxy - Fyz + C)*r(2) +...
        (Ez - Fxz - Fyz + C)*r(3);
R(3)  = (Fxy - C)*r(1)*r(2) + (Fxz - C)*r(1)*r(3) + (Fyz - C)*r(3)*r(2);
R(4)  = C*prod(r);

