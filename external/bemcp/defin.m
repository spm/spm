function defin

% defin :
%
% Save all the principal parameter in the file defin.mat
% 		n ndip_rad r1 r2 r3 sig1 sig2 sig3 sl_skull


fact=1 ; 	% 90 for more realistic head.

% Surface
%%%%%%%%%
n=32 % ATTENTION must be even !!!!!

% Npoints = 5*n^2/4+2 
% Ntri    = 5*n^2/2

% 8  -> 82   points -> 160 triangles
% 10 -> 127  points -> 250 triangles
% 12 -> 182  points -> 360 triangles
% 14 -> 247  points -> 490 triangles
% 16 -> 322  points -> 640 triangles
% 18 -> 407  points -> 810 triangles
% 20 -> 502  points -> 1000 triangles
% 22 -> 607  points -> 1210 triangles
% 24 -> 722  points -> 1440 triangles
% 26 -> 847  points -> 1690 triangles
% 28 -> 982  points -> 1960 triangles
% 30 -> 1127 points -> 2250 triangles
% 32 -> 1282 points -> 2560 triangles
% 34 -> 1447 points -> 2890 triangles
% 36 -> 1622 points -> 3240 triangles
% 40 -> 2002 points -> 4000 triangles
% 42 -> 2207 points -> 4410 triangles
% 44 -> 2422 points -> 4840 triangles
% 46 -> 2647 points -> 5290 triangles
% 48 -> 2882 points -> 5760 triangles
% 54 -> 3647 points -> 7290 triangles
% 56 -> 3922 points -> 7840 triangles
% 58 -> 4207 points -> 8410 triangles
% 60 -> 4502 points -> 9000 triangles
% 62 -> 4807 points -> 9610 triangles
% 64 -> 5122 points -> 10240 triangles.


% Volume
%%%%%%%%
ndip_rad = 10 % even , no more than 52 !!!

% 4  ->   142 dipole locations
% 6  ->   565 dipole locations
% 8  ->  1448 dipole locations
% 10 ->  2983 dipole locations
% 12 ->  5262 dipole locations
% 14 ->  8593 dipole locations
% 16 -> 13048 dipole locations
% 18 -> 18743 dipole locations
% 20 -> 25958 dipole locations
% 22 -> 34905 dipole locations
% 24 -> 45484 dipole locations
% 26 -> 58339 dipole locations
% 28 -> 73118 dipole locations
% 30 -> 90185 dipole locations


r1 = .8*fact ; r2 = .9*fact ; r3 = 1.*fact ;

%r1 = .9*fact ; r2 = .95*fact ; r3 = 1.*fact ;
%r1 = .6*fact ; r2 = .8*fact ; r3 = 1.*fact ;
%r1 = r3 ; r2 = r3 ;
%r1 = .4*fact ; r2 = .7*fact ; r3 = 1.*fact ;

% sig1 = .33 ; sig2 = .0042 ; sig3 = .33 ;
sig1 = 1 ; sig2 = .01 ; sig3 = 1 ;
% sig1 = .5 ; sig2 = 1 ; sig3 = 2 ;

sl_skull = 55 ; % inferior slice of usefull skull image 

cd /home/chrisp/matlab/Trian
save defin n ndip_rad r1 r2 r3 sig1 sig2 sig3 sl_skull

return

