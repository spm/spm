%function prep_C(electr,head,dipoles)
%
% function prep_C(electr,head,dipoles)
%
% Prepare the matrices :
%	- splitting and reordering lists
%	- weights
% Calculates the matrices Cij
%
% Save the results in param.mat, Cp1, Cp2
%
% Use calc_IFS then calc_L afterwards
%
% Last modified by chrisp@fil, 99/03/12 : Deflation for C.3 only !


load defin

tsph1 = head(1)
tsph2 = head(2)
tsph3 = head(3)
dsph  = dipoles

nvert = [tsph1.nr(1) tsph2.nr(1) tsph3.nr(1)] 
Max_nvert = max([tsph1.nr(1) tsph2.nr(1) tsph3.nr(1)]) 
nvert(find(nvert==Max_nvert)) = [] ;
Ma2 = max(nvert)

n_shareG = ceil(3*dsph.nr(1)/Ma2) 
shareG = ceil(dsph.nr(1)/n_shareG) ;
for i=1:n_shareG
	if i~=n_shareG
		eval(['share', num2str(i) ,'=(i-1)*shareG+1:i*shareG ;']) ;
	else
		eval(['share', num2str(i) ,'=(i-1)*shareG+1:dsph.nr(1) ;']) ;
	end
	eval(['lshare', num2str(i) ,'=length(share', num2str(i) ,') ']) ;
end

% Deflation for C.3 only !
Ntot_vert = tsph1.nr(1) + tsph2.nr(1) + tsph3.nr(1)
defl = 1/tsph3.nr(1)

Nel  = electr.nr(1) ;
El_vert = tsph3.tri(electr.ctri,:)' ;
El_vert = El_vert(:) ;
nEl_vert = 3*Nel ; % or length(El_vert) ;
Rem_vert = (1:tsph3.nr(1))' ;
Rem_vert(El_vert) = [] ;
nRem_vert = tsph3.nr(1)-nEl_vert ; % or length(Rem_vert) ;
ordered_list = [Rem_vert ; El_vert] ;

Perm_r = sparse(tsph3.nr(1),tsph3.nr(1)) ;
for i=1:tsph3.nr(1)
	Perm_r(i,ordered_list(i)) = 1 ;
end
Perm_c = Perm_r' ;

% Initialize weights :
weight(1) = (sig1-sig2)/((sig1+sig2)*2*pi) ;
weight(4) = (sig2-sig3)/((sig1+sig2)*2*pi) ;
weight(7) = sig3/((sig1+sig2)*2*pi) ;
weight(2) = (sig1-sig2)/((sig2+sig3)*2*pi) ;
weight(5) = (sig2-sig3)/((sig2+sig3)*2*pi) ;
weight(8) = sig3/((sig2+sig3)*2*pi) ;
weight(3) = (sig1-sig2)/(sig3*2*pi) ;
weight(6) = (sig2-sig3)/(sig3*2*pi) ;
weight(9) = 1/(2*pi) ;
weight(10) = 1/(2*pi*(sig1+sig2)) ;
weight(11) = 1/(2*pi*(sig2+sig3)) ;
weight(12) = 1/(2*pi*sig3) ;


clear dipoles* fact n ndip_rad r1 r2 r3 sl_skull
save Param_dip5

% No deflation for these matrices !
disp('flag 2.2')
C11st = Cii_lin(tsph1.tri,tsph1.vert,weight(1),0,tsph1.pt4.vert) ;
C21 = Cij_lin(tsph2.vert,tsph1.vert,tsph1.tri,weight(2),0) ;
save C11_21 C11st C21
clear C11st C21

disp('flag 2.3')
C22st = Cii_lin(tsph2.tri,tsph2.vert,weight(5),0,tsph2.pt4.vert) ;
C12 = Cij_lin(tsph1.vert,tsph2.vert,tsph2.tri,weight(4),0) ;
save C12_22 C12 C22st
clear C12 C22st

disp('flag 2.4')
C31 = Perm_r*Cij_lin(tsph3.vert,tsph1.vert,tsph1.tri,weight(3),0) ;
save C31 C31
clear C31

C32 = Perm_r*Cij_lin(tsph3.vert,tsph2.vert,tsph2.tri,weight(6),0) ;
save C32 C32
clear C32

% Only deflation here !
disp('flag 2.5')
C13 = Cij_lin(tsph1.vert,tsph3.vert,tsph3.tri,weight(7),defl)*Perm_c ;
C23 = Cij_lin(tsph2.vert,tsph3.vert,tsph3.tri,weight(8),defl)*Perm_c ;
C33st = Perm_r * ...
	Cii_lin(tsph3.tri,tsph3.vert,weight(9),defl,tsph3.pt4.vert)*Perm_c ;

save C13_23_33 C13 C23 C33st
