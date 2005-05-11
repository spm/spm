function [el_sphc,el_name] = spm_eegfp_electrset(el_set)
%------------------------------------------------------------------------
% FORMAT [el_sphc,el_name] = spm_eegfp_electrset(el_set) ;
% or
% FORMAT [set_Nel,set_name] = spm_eegfp_electrset ;
% 
% Creates the electrode set on a sphere, set type is defined by 'el_set'.
% or return the name of the sets available, and the number of electrodes.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips,
% $Id$

set_name = strvcat('10-20 system with 23 electrodes.', ...
    '10-20 system with 19 electrodes.', ...
    '61 Equidistant electrodes.', ...
    '148 electrodes, Pascual-Marqui''s set.', ...
    '31 electrodes, Rik''s set.', ...
    '10-20 system with 21 electrodes.', ...
    '10-20 system with 29 electrodes.', ...
    '62 electrodes on ext. 10-20 system.');
Nr_electr = [23 19 61 148 31 21 29 62];
if nargin<1
    el_sphc = Nr_electr;
    el_name = set_name;
else
     try
		switch el_set
			case 1
				fprintf(['\n',set_name(el_set,:),'\n\n']);
				[el_sphc,el_name] = sys1020_23 ;
			case 2
				fprintf(['\n',set_name(el_set,:),'\n\n']);
				[el_sphc,el_name] = sys1020_19 ;
			case 3
				fprintf(['\n',set_name(el_set,:),'\n\n']);
				[el_sphc,el_name] = sys61 ;
			case 4
   				fprintf(['\n',set_name(el_set,:),'\n\n']);
				[el_sphc,el_name] = sys148 ;
			case 5
   				fprintf(['\n',set_name(el_set,:),'\n\n']);
				el_rik = [59  1 14 50 48 33 19 47 31 17 ...
                          46 30 29 44 36  9 22 38 11 24 ...
                          39 26 25 40 42 53 37 49 41 45  8];
                % Subset of electrodes for Rik's expe.
                [el_sphc,el_name] = sys61 ;
				el_sphc = el_sphc(:,el_rik) ;
                el_name = el_name(el_rik,:) ;
			case 6
				fprintf(['\n',set_name(el_set,:),'\n\n']);
                el_not_21 = [2 22];
                el_21 = 1:23; el_21(el_not_21) = [];
                [el_sphc,el_name] = sys1020_23 ;
				el_sphc = el_sphc(:,el_21) ;
                el_name = el_name(el_21,:) ;
			case 7
				fprintf(['\n',set_name(el_set,:),'\n\n']);
				[el_sphc,el_name] = sys1020_29 ;
            case 8
				fprintf(['\n',set_name(el_set,:),'\n\n']);
				[el_sphc,el_name] = sys1020_62 ;               
			otherwise
				warning('Unknown electrodes set !!!')
                el_sphc = []; el_name = [];
		end
        el_name = cellstr(el_name);
    catch
        error(['Provide the electrode set number : from 1 to ',num2str(length(Nr_electr))])
    end
end

%________________________________________________________________________
%
% Subfunctions defining the electrode sets
%
%________________________________________________________________________
% FORMAT el_sphc = sysXXXX ;
%
% Genrates a set of electrode coordinates on a sphere.
% Other sets can be added here, as long as the format is respected :
% 	el_sphc	: 3xNel coordinates on a radius 1 sphere
%	Nel	: nr of electrodes
%------------------------------------------------------------------------
function [el_sphc,el_name] = sys1020_23 ;

% Classic 10-20 system with 23 usefull electrodes
% The mastoid electrodes (A1-A2) are located on the hemispheric plane

d2r=pi/180; % converts degrees to radian
Fpz = [ 0 sin(72*d2r) cos(72*d2r) ] ;
Fp1 = [ -sin(72*d2r)*sin(18*d2r) sin(72*d2r)*cos(18*d2r) cos(72*d2r) ] ;
F7 = [ -sin(72*d2r)*sin(54*d2r) sin(72*d2r)*cos(54*d2r) cos(72*d2r) ] ;
T3 = [ -sin(72*d2r) 0 cos(72*d2r) ] ;
C3 = [ -sin(36*d2r) 0 cos(36*d2r) ] ;
Cz = [ 0 0 1 ] ;
Fz = [ 0 sin(36*d2r) cos(36*d2r) ] ;
temp=Fz+F7 ;
F3=temp/norm(temp);

el_sphc=[ Fp1 ;						% Fp1
	Fpz						% Fpz
	-Fp1(1)	 	Fp1(2) 		Fp1(3) ;	% Fp2
	F7 ;						% F7
	F3 ;						% F3
	Fz ;						% Fz
	-F3(1)		F3(2)		F3(3) ;		% F4
	-F7(1)		F7(2) 		F7(3) ;		% F8
	-1		0		0 ;		% A1
	T3 ;						% T3
	C3 ;						% C3
	Cz ;						% Cz
	-C3(1)		C3(2)		C3(3) ;		% C4
	-T3(1)		T3(2)		T3(3) ;		% T4
	1		0		0 ;		% A2
	F7(1)		-F7(2)		F7(3) ;		% T5
	F3(1)		-F3(2)		F3(3) ;		% P3
	Fz(1)		-Fz(2) 		Fz(3) ;		% Pz
	-F3(1)		-F3(2)		F3(3) ;		% P4
	-F7(1)		-F7(2)		F7(3) ;		% T6
	Fp1(1)	 	-Fp1(2)	 	Fp1(3) ;	% O1
	Fpz(1)		-Fpz(2) 	Fpz(3) ;	% Oz
	-Fp1(1)	 	-Fp1(2) 	Fp1(3) ]' ;	% O2
el_name = str2mat('Fp1','Fpz','Fp2', ...
		  'F7','F3','Fz','F4','F8', ...
		  'A1','T3','C3','Cz','C4','T4','A2', ...
		  'T5','P3','Pz','P4','T6', ...
		  'O1','Oz','O2') ;
return
%------------------------------------------------------------------------
function [el_sphc,el_name] = sys1020_19 ;

% Classic 10-20 system with 19 usefull electrodes
% The mastoid electrodes (A1-A2) would be located on the hemispheric plane
% Compared with the previous set, electrodes Fpz, A1, A2 & Oz are removed

d2r=pi/180; % converts degrees to radian
Fpz = [ 0 sin(72*d2r) cos(72*d2r) ] ;
Fp1 = [ -sin(72*d2r)*sin(18*d2r) sin(72*d2r)*cos(18*d2r) cos(72*d2r) ] ;
F7 = [ -sin(72*d2r)*sin(54*d2r) sin(72*d2r)*cos(54*d2r) cos(72*d2r) ] ;
T3 = [ -sin(72*d2r) 0 cos(72*d2r) ] ;
C3 = [ -sin(36*d2r) 0 cos(36*d2r) ] ;
Cz = [ 0 0 1 ] ;
Fz = [ 0 sin(36*d2r) cos(36*d2r) ] ;
temp=Fz+F7 ;
F3=temp/norm(temp);

el_sphc=[ Fp1 ;						% Fp1
	-Fp1(1)	 	Fp1(2) 		Fp1(3) ;	% Fp2
	F7 ;						% F7
	F3 ;						% F3
	Fz ;						% Fz
	-F3(1)		F3(2)		F3(3) ;		% F4
	-F7(1)		F7(2) 		F7(3) ;		% F8
	T3 ;						% T3
	C3 ;						% C3
	Cz ;						% Cz
	-C3(1)		C3(2)		C3(3) ;		% C4
	-T3(1)		T3(2)		T3(3) ;		% T4
	F7(1)		-F7(2)		F7(3) ;		% T5
	F3(1)		-F3(2)		F3(3) ;		% P3
	Fz(1)		-Fz(2) 		Fz(3) ;		% Pz
	-F3(1)		-F3(2)		F3(3) ;		% P4
	-F7(1)		-F7(2)		F7(3) ;		% T6
	Fp1(1)	 	-Fp1(2)	 	Fp1(3) ;	% O1
	-Fp1(1)	 	-Fp1(2) 	Fp1(3) ]' ;	% O2
el_name = str2mat('Fp1','Fp2', ...
		  'F7','F3','Fz','F4','F8', ...
		  'T3','C3','Cz','C4','T4', ...
		  'T5','P3','Pz','P4','T6', ...
		  'O1','O2') ;
return
%------------------------------------------------------------------------
function [el_sphc,el_name] = sys1020_29 ;

% Classic 10-20 system with 29 useful electrodes
% The mastoid electrodes (A1-A2) are located on the hemispheric plane
% Intermediate locations at 10% are used here, hence a few "classic"
% electrodes have different names (EG. T8=T4, P8=T6, T7=T3, P7=T5

d2r=pi/180; % converts degrees to radian
Fpz = [ 0 sin(72*d2r) cos(72*d2r) ] ;
Fp1 = [ -sin(72*d2r)*sin(18*d2r) sin(72*d2r)*cos(18*d2r) cos(72*d2r) ] ;
F7 = [ -sin(72*d2r)*sin(54*d2r) sin(72*d2r)*cos(54*d2r) cos(72*d2r) ] ;
T3 = [ -sin(72*d2r) 0 cos(72*d2r) ] ;
C3 = [ -sin(36*d2r) 0 cos(36*d2r) ] ;
Cz = [ 0 0 1 ] ;
Fz = [ 0 sin(36*d2r) cos(36*d2r) ] ;
F3 = Fz+F7 ;F3 = F3/norm(F3);

Fc5 = F7+F3+T3+C3; Fc5 = Fc5/norm(Fc5);
Fc1 = F3+Fz+C3+Cz; Fc1 = Fc1/norm(Fc1);

el_sphc=[ 	1		0		0 ;		    % A2
	-F7(1)		F7(2) 		F7(3) ;		% F8
	-T3(1)		T3(2)		T3(3) ;		% T8=T4 or flipped T3
	-F7(1)		-F7(2)		F7(3) ;		% P8=T6
    -Fc5(1)     Fc5(2)      Fc5(3) ;    % Fc6
    -Fc5(1)     -Fc5(2)     Fc5(3) ;    % Cp6
	-Fp1(1)	 	Fp1(2) 		Fp1(3) ;	% Fp2
	-F3(1)		F3(2)		F3(3) ;		% F4
	-C3(1)		C3(2)		C3(3) ;		% C4
   	-F3(1)		-F3(2)		F3(3) ;		% P4
 	-Fp1(1)	 	-Fp1(2) 	Fp1(3) ;	% O2
    -Fc1(1)     Fc1(2)      Fc1(3) ;    % Fc2
    -Fc1(1)     -Fc1(2)      Fc1(3) ;   % Cp2
	Fz ;						        % Fz
	Cz ;						        % Cz
	Fz(1)		-Fz(2) 		Fz(3) ;		% Pz
    Fc1                                 % Fc1
    Fc1(1)     -Fc1(2)      Fc1(3) ;    % Cp1
    Fp1 ;						        % Fp1
    F3 ;						        % F3
	C3 ;						        % C3
	F3(1)		-F3(2)		F3(3) ;		% P3
	Fp1(1)	 	-Fp1(2)	 	Fp1(3) ;	% O1
    Fc5 ;                               % Fc5
    Fc5(1)     -Fc5(2)     Fc5(3) ;     % Cp5
	F7 ;						        % F7
	T3 ;						        % T7=T3
	F7(1)		-F7(2)		F7(3) ;		% P7=T5
	-1		    0		    0 ]' ;      % A1

el_name = strvcat('A2','F8','T8','P8','Fc6','Cp6', ...
          'Fp2','F4','C4','P4','O2','Fc2','Cp2', ...
          'Fz','Cz','Pz','Fc1','Cp1', ...
          'Fp1','F3','C3','P3','O1','Fc5','Cp5', ...
          'F7','T7','P7','A1') ;
return
%------------------------------------------------------------------------
function [el_sphc,el_name] = sys61 ;

% 61  quasi-equidistant electrodes, as used on the easycap.
%

d2r = pi/180 ;
Ne_pr = [ 6 12 15 16 14 ] ;
	% Number of electrodes per "ring", from top to bottom
	% Cz is assumed to be perfectly on top [0 0 1]

N_r = length(Ne_pr) ;
dth = 90 * d2r / N_r ;

el_sphc = [0 0 1] ;
for i=1:N_r
	thet = i * dth ;
	sth = sin(thet) ;
	cth = cos(thet) ;
	l_phi = (90:-360/Ne_pr(i):-269)*d2r ;
	for phi=l_phi
		el_sphc = [ el_sphc ; sth*cos(phi) sth*sin(phi) cth ] ;
	end
end

% remove 3 frontal electrodes at the level of eyes.
el_sphc([51 52 64],:) = [] ;
el_sphc = el_sphc';
el_name = num2str((1:61)');
%figure,plot3(el_sphc(1,:),el_sphc(2,:),el_sphc(3,:),'*');axis vis3d;rotate3d on
return
%------------------------------------------------------------------------
function [el_sphc,el_name] = sys148 ;

% 148 electrodes spread as in the paper by Pascual-Marqui et al.
% This si not suitable for a realist head model, as the location of the eyes
% is ignored...

warning('This electrode set is not suitable for realistic head model !');

d2r = pi/180 ;
Ne_pr = [ 6 12 17 21 23 24 23 21 ] ;
el_sphc = [0 0 1] ;
N_r = length(Ne_pr) ;
dth = 90 * d2r / (1+length(find(diff(Ne_pr)>0))) ;
	% Number of electrodes per "ring", from top to bottom
	% The nr of electr per ring decreases under the hemispheric line
	% Cz is assumed to be perfectly on top [0 0 1]

for i=1:N_r
	thet = i * dth ;
	sth = sin(thet) ;
	cth = cos(thet) ;
	dphi = - 360 / Ne_pr(i) * d2r ;
	phi0 = (1-rem(i,2)) * dphi / 2 ;
	for j=1:Ne_pr(i)
		phi = phi0 + (j-1)*dphi ;
		el_sphc = [ el_sphc ; -sth*sin(phi) sth*cos(phi) cth ] ;
	end
end
el_sphc = el_sphc';
el_name = num2str((1:148)');
%figure,plot3(el_sphc(1,:),el_sphc(2,:),el_sphc(3,:),'*');axis vis3d;rotate3d on
return

%------------------------------------------------------------------------
function [el_sphc,el_name] = sys1020_62 ;

% Classic 10-20 system with 62 useful electrodes
% The mastoid electrodes (A1-A2) are located on the hemispheric plane
% Intermediate locations at 10% are used here, hence a few "classic"
% electrodes have different names (eg. T8=T4, P8=T6, T7=T3, P7=T5

d2r=pi/180; % converts degrees to radian

% Electrodes on the front left quarter
%-------------------------------------
% Main electrodes
FPz = [ 0 sin(72*d2r) cos(72*d2r) ] ;
FP1 = [ -sin(72*d2r)*sin(18*d2r) sin(72*d2r)*cos(18*d2r) cos(72*d2r) ] ;
AF7 = [ -sin(72*d2r)*sin(36*d2r) sin(72*d2r)*cos(36*d2r) cos(72*d2r) ] ;
F7  = [ -sin(72*d2r)*sin(54*d2r) sin(72*d2r)*cos(54*d2r) cos(72*d2r) ] ;
FT7 = [ -sin(72*d2r)*sin(72*d2r) sin(72*d2r)*cos(72*d2r) cos(72*d2r) ] ;
T7  = [ -sin(72*d2r) 0 cos(72*d2r) ] ;
C5  = [ -sin(54*d2r) 0 cos(54*d2r) ] ;
C3  = [ -sin(36*d2r) 0 cos(36*d2r) ] ;
C1  = [ -sin(18*d2r) 0 cos(18*d2r) ] ;
AFz = [ 0 sin(54*d2r) cos(54*d2r) ] ;
Fz  = [ 0 sin(36*d2r) cos(36*d2r) ] ;
Ref = [ 0 sin(18*d2r) cos(18*d2r) ] ;
Cz  = [ 0 0 1 ] ;
TP9 = [ -sin(108*d2r) cos(108*d2r) 0 ] ;

% Intermediate electrodes
AF3 = AF7+AFz ; AF3 = AF3/norm(AF3);
F3  = Fz+F7   ; F3  = F3 /norm(F3);
FC3 = FT7+Ref ; FC3 = FC3/norm(FC3);
F5  = F7+F3   ; F5  = F5 /norm(F5);
F1  = Fz+F3   ; F1  = F1 /norm(F1);
FC5 = FT7+FC3 ; FC5 = FC5/norm(FC5);
FC1 = Ref+FC3 ; FC1 = FC1/norm(FC1);

% Note, there is no A1/A2 inthis montage but TP9 and TP10

el_sphc=[ FP1 ;                         % Fp1
    FPz ;                               % Fpz
	-FP1(1)	 	FP1(2) 		FP1(3) ;	% Fp2
    AF7 ;                               % AF7
    AF3 ;                               % AF3
    AFz ;                               % AFz
    -AF3(1)     AF3(2)      AF3(3) ;    % AF4
    -AF7(1)     AF7(2)      AF7(3) ;    % AF8
    F7 ;                                % F7
    F5 ;                                % F5
    F3 ;						        % F3
    F1 ;						        % F1
    Fz ;						        % Fz
	-F1(1)	 	F1(2) 		F1(3) ;	    % F2
    -F3(1)      F3(2)       F3(3) ;     % F4
    -F5(1)      F5(2)       F5(3) ;     % F6
    -F7(1)      F7(2)       F7(3) ;     % F8
    FT7 ;                               % FT7
    FC5 ;                               % FC5
    FC3 ;                               % FC3
    FC1 ;                               % FC1
	-FC1(1)	 	FC1(2) 		FC1(3) ;	% FC2
	-FC3(1)	 	FC3(2) 		FC3(3) ;	% FC4
	-FC5(1)	 	FC5(2) 		FC5(3) ;	% FC6
	-FT7(1)	 	FT7(2) 		FT7(3) ;	% FT8
    T7 ;                                % T7
    C5 ;                                % C5
    C3 ;                                % C3
    C1 ;                                % C1
	Cz ;						        % Cz
	-C1(1)		C1(2)		C1(3) ;		% C2
	-C3(1)		C3(2)		C3(3) ;		% C4
	-C5(1)		C5(2)		C5(3) ;		% C6
    -T7(1)      T7(2)       T7(3) ;     % T8    Start posteririor part
    TP9 ;                               % TP9
	FT7(1)	 	-FT7(2)     FT7(3) ;	% TP7
	FC5(1)	 	-FC5(2)		FC5(3) ;	% CP5
	FC3(1)	 	-FC3(2) 	FC3(3) ;	% CP3
	FC1(1)	 	-FC1(2)		FC1(3) ;	% CP1
    Ref(1)      -Ref(2)     Ref(3) ;    % Cpz
	-FC1(1)	 	-FC1(2)		FC1(3) ;	% CP2
	-FC3(1)	 	-FC3(2) 	FC3(3) ;	% CP3
	-FC5(1)	 	-FC5(2)		FC5(3) ;	% CP6
	-FT7(1)	 	-FT7(2)     FT7(3) ;	% TP8
    -TP9(1)     TP9(2)      TP9(3) ;    % TP10
    F7(1)       -F7(2)       F7(3) ;    % P7
    F5(1)       -F5(2)       F5(3) ;    % P5
    F3(1)       -F3(2)       F3(3) ;    % P3
	F1(1)	 	-F1(2) 		F1(3) ;	    % P1
	Fz(1)		-Fz(2) 		Fz(3) ;		% Pz 
 	-F1(1)	 	-F1(2) 		F1(3) ;	    % P2
    -F3(1)      -F3(2)       F3(3) ;    % P4
    -F5(1)      -F5(2)       F5(3) ;    % P6
    -F7(1)      -F7(2)       F7(3) ;    % P8
    AF7(1)      -AF7(2)     AF7(3) ;    % PO7
    AF3(1)      -AF3(2)     AF3(3) ;    % PO3
    AFz(1)      -AFz(2)     AFz(3) ;    % POz
    -AF3(1)     -AF3(2)     AF3(3) ;    % PO4
    -AF7(1)      -AF7(2)     AF7(3) ;    % PO8
	FP1(1)	 	-FP1(2)	 	FP1(3) ;	% O1
	FPz(1)	 	-FPz(2)	 	FPz(3) ;	% Oz
 	-FP1(1)	 	-FP1(2) 	FP1(3) ]';	% O2

el_name = strvcat('FP1','FPz','FP2','AF7','AF3','AFz','AF4','AF8', ...
                  'F7','F5','F3','F1','Fz','F2','F4','F6','F8', ...
                  'FT7','FC5','FC3','FC1','FC2','FC4','FC6','FT8', ...
                  'T7','C5','C3','C1','Cz','C2','C4','C6','T8', ...
                  'TP9','TP7','CP5','CP3','CP1','CPz','CP2','CP4','CP6','TP8','TP10', ...
                  'P7','P5','P3','P1','Pz','P2','P4','P6','P8', ...
                  'PO7','PO3','POz','PO4','PO8','O1','Oz','O2');

return
%------------------------------------------------------------------------
