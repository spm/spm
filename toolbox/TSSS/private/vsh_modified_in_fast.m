% Vm = vsh_modified_in_fast(theta,phi,l,m,p0)
%
% Vector spherical harmonic funtion modified from Hill's function
% V_lm for quantum numbers l and m. Angles theta and phi should be
% given in radians. p0 contains the values of the Legendre polynomial.
%
function Vm = vsh_modified_in_fast(theta,phi,l,m,p0)

% Copyright (c) 2016, Elekta Oy
% ---------------------------------------
% 
% Redistribution and use of the Software in source and binary forms, with or without 
% modification, are permitted for non-commercial use.
% 
% The Software is provided "as is" without warranties of any kind, either express or
% implied including, without limitation, warranties that the Software is free of defects,
% merchantable, fit for a particular purpose. Developer/user agrees to bear the entire risk 
% in connection with its use and distribution of any and all parts of the Software under this license.
% 

scale = 1;

%scale_sph = (-1)^m*sqrt((2*l+1)*prod(1:(l-m))/(4*pi*prod(1:(l+m))));
%scale_minus = 1/((-1)^(m-1)*sqrt((2*l+1)*prod(1:(l-(m-1)))/(4*pi*prod(1:(l+(m-1))))));
%scale_plus = 1/((-1)^(m+1)*sqrt((2*l+1)*prod(1:(l-(m+1)))/(4*pi*prod(1:(l+(m+1))))));

Y = spharm_fast(theta,phi,l,m,p0);
if m > -l
  Yminus = spharm_fast(theta,phi,l,m-1,p0);
else
  Yminus = 0;
end
if m < l
  Yplus = spharm_fast(theta,phi,l,m+1,p0);
else
  Yplus = 0;
end

cosp = cos(phi); sinp = sin(phi);
%dY = 0.5*scale_sph*((l+m)*(l-m+1)*scale_minus*Yminus*complex(cos(phi),sin(phi))-scale_plus*Yplus*complex(cos(-phi),sin(-phi)));  % dY/dtheta
%dY = 0.5*(-sqrt((l+m)*(l-m+1))*Yminus*complex(cos(phi),sin(phi))+sqrt((l-m)*(l+m+1))*Yplus*complex(cos(-phi),sin(-phi)));
dY = 0.5*(-sqrt((l+m)*(l-m+1))*Yminus*complex(cosp,sinp)+sqrt((l-m)*(l+m+1))*Yplus*complex(cosp,-sinp));
if theta == 0 % To prevent division by zero
  theta = eps;
end
Vm(1) = scale*(-(l+1))*Y;
Vm(2) = scale*dY;
Vm(3) = scale*(i*m/sin(theta))*Y;



