%
% [Sout,SNout] = Sout_vsh_fast(r_sphere,R,EX,EY,EZ,ch_types,Lout)
%
% Calculate the outer SSS basis Sout using
% vector spherical harmonics
%
function [Sout,SNout] = Sout_vsh_fast(r_sphere,R,EX,EY,EZ,ch_types,Lout)

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

MAG = 1;
GRAD = 0;
mu0 = 1.25664e-6; % Permeability of vacuum
%
% For numerical surface integration:
%
mag_size = 21e-3;
baseline = 16.69e-3;
%d = sqrt(3/5)*mag_size/2;
%dx1 = 5.89e-3;
%dx2 = 10.8e-3;
%dy = 6.71e-3;
%Dmag = [0 0; d d; -d d; -d -d; d -d; 0 d; 0 -d; d 0; -d 0]';
%Dgrad = [dx1 dy; dx2 dy; dx1 -dy; dx2 -dy; -dx1 dy; -dx2 dy; -dx1 -dy; -dx2 -dy]';
%weights_mag = [16/81 25/324 25/324 25/324 25/324 10/81 10/81 10/81 10/81]';
%for j = 1:8
%  if j <= 4
%    weights_grad(j) = 1/(4*baseline);
%  else
%    weights_grad(j) = -1/(4*baseline);
%  end
%end
%weights_grad = weights_grad';
%
% 27.11.09 STa: The following integration points correspond to the normal coil
% accuracy in Xfit
%
Dmag = 1e-3*[6.45 6.45; 6.45 -6.45; -6.45 6.45; -6.45 -6.45]';
weights_mag = (1/4)*ones(4,1);
Dgrad = 1e-3*[8.4 0; -8.4 0]';
weights_grad = (1/baseline)*[1;-1];

nchan = length(ch_types);
for ch = 1:nchan
  %disp(ch)
  count = 1;
  R(:,ch) = R(:,ch) - r_sphere;
  if ch_types(ch) == GRAD
    D = Dgrad;
    weights = weights_grad;
  elseif ch_types(ch) == MAG
    D = Dmag;
    weights = weights_mag;
  else
    error('Unknown sensor type!');
  end
  Sout(ch,:) = -mu0*vsh_response(R(:,ch),EX(:,ch),EY(:,ch),EZ(:,ch),D,weights,Lout);
  %for l = 1:Lout
  %  for m = -l:l
  % Sout(ch,count) = -mu0*vsh_response(R(:,ch),EX(:,ch),EY(:,ch),EZ(:,ch),D,weights,l,m);
  %count = count + 1;
  %end
  %end
end
for j = 1:size(Sout,2)
  SNout(:,j) = Sout(:,j)/norm(Sout(:,j));
end


function Sout_elements = vsh_response(r,ex,ey,ez,D,weights,Lout)

for j = 1:length(weights)
  r_this = r + D(1,j)*ex + D(2,j)*ey;
  rn(j) = norm(r_this);
  theta(j) = acos(r_this(3)/rn(j));
  phi(j) = atan2(r_this(2),r_this(1));
  sint(j) = sin(theta(j));
  sinp(j) = sin(phi(j));
  cost(j) = cos(theta(j));
  cosp(j) = cos(phi(j));
  for l = 1:Lout
    p0{l}(:,j) = legendre(l,cos(theta(j)));
    rnv(j,l) = rn(j)^(l-1);
  end
end
Sout_elements = [];
for l = 1:Lout
  for m = -l:l
    for j = 1:length(weights)
      %disp(j)
      ws = vsh_modified_out_fast(theta(j),phi(j),l,m,p0{l}(:,j))'*(rnv(j,l));
      W(1,j) = ws(1)*sint(j)*cosp(j) + ws(2)*cost(j)*cosp(j) - ws(3)*sinp(j);
      W(2,j) = ws(1)*sint(j)*sinp(j) + ws(2)*cost(j)*sinp(j) + ws(3)*cosp(j);
      W(3,j) = ws(1)*cost(j) - ws(2)*sint(j);
    end
    Sout_elements = [Sout_elements dot(W*weights,ez)]; % Cartesian coordinates
                                                     %Sin_element = Sin_element/sqrt((l+1)*(2*l+1));  % Back to orthonormal presentation
  end 
end 




