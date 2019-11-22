%
% [Sout,SNout] = Sout_vsh_ctf_fast(r_sphere,R,EX,EY,EZ,Lout)
%
% Calculate the external SSS basis Sout for a CTF system
% using vector spherical harmonics
%
function [Sout,SNout] = Sout_vsh_ctf_fast(r_sphere,R,EX,EY,EZ,Lout)

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
   count = 1;
   R(:,ch) = R(:,ch) - r_sphere;
   Sout(ch,:) = -mu0*vsh_response(R(:,ch),EX(:,ch),EY(:,ch),EZ(:,ch),D,weights,Lout);
end
for j = 1:size(Sout,2)
   SNout(:,j) = Sout(:,j)/norm(Sout(:,j));
end

Sout  = Sout(:, 4:end);
SNout = SNout(:, 4:end);


function Sout_elements = vsh_response(r,ex,ey,ez,D,weights,Lin)

for j = 1:length(weights)
    r_this = r + D(1,j)*ex + D(2,j)*ey;
    rn(j) = norm(r_this);
    theta(j) = acos(r_this(3)/rn(j));
    phi(j) = atan2(r_this(2),r_this(1));
    sint(j) = sin(theta(j));
    sinp(j) = sin(phi(j));
    cost(j) = cos(theta(j));
    cosp(j) = cos(phi(j));
    for l = 1:Lin
       p0{l}(:,j) = legendre(l,cos(theta(j)));
       rnv(j,l) = rn(j)^(l+2);
   end
end
Sout_elements = [];
for l = 1:Lin
  for m = -l:l
    for j = 1:length(weights)
      %vs = vsh_modified_in_fast(theta(j),phi(j),l,m,p0{l}(:,j))'/rn(j)^(l+2);
      vs = vsh_modified_out_fast(theta(j),phi(j),l,m,p0{l}(:,j))'/rnv(j,l);
      V(1,j) = vs(1)*sint(j)*cosp(j) + vs(2)*cost(j)*cosp(j) - vs(3)*sinp(j);
      V(2,j) = vs(1)*sint(j)*sinp(j) + vs(2)*cost(j)*sinp(j) + vs(3)*cosp(j);
      V(3,j) = vs(1)*cost(j) - vs(2)*sint(j);
    end
    Sout_elements = [Sout_elements dot(V*weights,ez)]; % Cartesian coordinates
                                                     %Sin_element = Sin_element/sqrt((l+1)*(2*l+1));  % Back to orthonormal presentation
  end 
end 
%Sout_element = dot(V*weights,ez); % Cartesian coordinates



