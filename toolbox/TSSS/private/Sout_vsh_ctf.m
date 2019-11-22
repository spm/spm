%
% [Sout,SNout] = Sout_vsh_ctf(r_sphere,R,EX,EY,EZ,Lout)
%
% Calculate the external SSS basis Sout for a CTF system
% using vector spherical harmonics
%
function [Sout,SNout] = Sout_vsh_ctf(r_sphere,R,EX,EY,EZ,Lin)

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

mu0 = 1.25664e-6; % Permeability of vacuum
%
% For numerical surface integration:
%
%baseline = 50e-3;
dx = 4.5e-3;
dy = 4.5e-3;
dz1 = 0;
dz2 = 50e-3;
D = [dx dy dz1; dx -dy dz1; -dx dy dz1; -dx -dy dz1; dx dy dz2; dx -dy dz2; -dx dy dz2; -dx -dy dz2]';
for j = 1:8
   if j <= 4
      weights(j) = 1/(4*1);
   else
      weights(j) = -1/(4*1);
   end
end
weights = weights';

for ch = 1:size(R,2)
   disp(ch)
   count = 1;
   R(:,ch) = R(:,ch) - r_sphere;
   for l = 1:Lin
      for m = -l:l
	 Sout(ch,count) = -mu0*vsh_response(R(:,ch),EX(:,ch),EY(:,ch),EZ(:,ch),D,weights,l,m);
	 count = count + 1;
      end
   end
end
for j = 1:size(Sout,2)
   SNout(:,j) = Sout(:,j)/norm(Sout(:,j));
end

Sout  = Sout(:, 4:end);
SNout = SNout(:, 4:end);

function Sout_element = vsh_response(r,ex,ey,ez,D,weights,l,m)

for j = 1:length(weights)
  r_this = r + D(1,j)*ex + D(2,j)*ey + D(3,j)*ez;
  rn = norm(r_this);
  theta = acos(r_this(3)/rn);
  phi = atan2(r_this(2),r_this(1));
  sint = sin(theta);
  sinp = sin(phi);
  cost = cos(theta);
  cosp = cos(phi);
  vs = vsh_modified_out(theta,phi,l,m)'*rn^(l-1);
  V(1,j) = vs(1)*sint*cosp + vs(2)*cost*cosp - vs(3)*sinp;
  V(2,j) = vs(1)*sint*sinp + vs(2)*cost*sinp + vs(3)*cosp;
  V(3,j) = vs(1)*cost - vs(2)*sint;
end
Sout_element = dot(V*weights,ez); % Cartesian coordinates



