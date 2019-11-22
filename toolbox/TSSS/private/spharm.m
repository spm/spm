% Y = spharm(theta,phi,l,m)
%
% Orthonormal spherical harmonic function (Arfken, p. 681) for quantum
% numbers l and m. Angles theta and phi should be given in radians.
%
function Y = spharm(theta,phi,l,m)

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

p0 = legendre(l,cos(theta));
p = p0(2:end);  % Values m > 1
%
% Arfken, p. 681. The Condon-Shortley phase is included in the
% Legendre polynomials
%
scale = sqrt((2*l+1)*prod(1:(l-m))/(4*pi*prod(1:(l+m)))); 
if m < 0
   Y = scale*((-1)^m)*(prod(1:(l-abs(m)))/prod(1:(l+abs(m))))*p(-m);
   ((-1)^m)*(prod(1:(l-m))/prod(1:(l+m)));
elseif m > 0
   Y = scale*p(m);
else
   Y = scale*p0(1);
end
Y = Y*complex(cos(m*phi),sin(m*phi));


