function L = bem_full_cst(head_model,dipoles)

% L = bem_full_cst([head_model,dipoles])
%
% Calculate the lead field for all the vertices of the 3rd surface
% uses the standard Cii, Cij, Gi function with the cst approximation.

TICTOC = clock ;

% 1. INITIALISATION
%%%%%%%%%%%%%%%%%%%
disp('flag 1')

load defin

if nargin==0
	n1 = n ; n2 = n ; n3 = n ;
%	n1 = 2*round(.75*n) ; n2 = n ; n3 = n ;

	% Sphere1 = brain
	tsph1 = tes_sph(n1,r1) ;
	tsph1 = m_edges(tsph1) ;
	tsph1 = tri_4thpt(tsph1)
	% Sphere2 = skull
	tsph2 = tes_sph(n2,r2) ;
	tsph2 = m_edges(tsph2) ;
	tsph2 = tri_4thpt(tsph2)
	% Sphere3 = scalp
	tsph3 = tes_sph(n3,r3) ;
	tsph3 = m_edges(tsph3) ;
	tsph3 = tri_4thpt(tsph3)
	% dipoles
	dsph=dip_sph(ndip_rad,r1)
else
	% With real head model !

end
Ntot_vert = tsph1.nr(1) + tsph2.nr(1) + tsph3.nr(1)
%defl = 1/Ntot_vert
%defl1 = defl ; defl2 = defl ; defl3 = defl ;
defl1 = 0 ; defl2 = 0 ; defl3 = 1/tsph3.nr(1) ;

% 2. MAIN MATRICES
%%%%%%%%%%%%%%%%%%

disp('flag 2.2')
weight = (sig1-sig2)/((sig1+sig2)*2*pi) ;
C11st = bem_Cii_cst(tsph1.tri,tsph1.vert,weight,defl1,tsph1.pt4.vert) ;
weight = (sig2-sig3)/((sig1+sig2)*2*pi) ;
C12 = bem_Cij_cst(tsph1.vert,tsph2.vert,tsph2.tri,weight,defl2) ;
weight = sig3/((sig1+sig2)*2*pi) ;
C13 = bem_Cij_cst(tsph1.vert,tsph3.vert,tsph3.tri,weight,defl3) ;

disp('flag 2.3')
weight = (sig1-sig2)/((sig2+sig3)*2*pi) ;
C21 = bem_Cij_cst(tsph2.vert,tsph1.vert,tsph1.tri,weight,defl1) ;
weight = (sig2-sig3)/((sig2+sig3)*2*pi) ;
C22st = bem_Cii_cst(tsph2.tri,tsph2.vert,weight,defl2,tsph2.pt4.vert) ;
weight = sig3/((sig2+sig3)*2*pi) ;
C23 = bem_Cij_cst(tsph2.vert,tsph3.vert,tsph3.tri,weight,defl3) ;

disp('flag 2.4')
weight = (sig1-sig2)/(sig3*2*pi) ;
C31 = bem_Cij_cst(tsph3.vert,tsph1.vert,tsph1.tri,weight,defl1) ;
weight = (sig2-sig3)/(sig3*2*pi) ;
C32 = bem_Cij_cst(tsph3.vert,tsph2.vert,tsph2.tri,weight,defl2) ;
weight = 1/(2*pi) ;
C33st = bem_Cii_cst(tsph3.tri,tsph3.vert,weight,defl3,tsph3.pt4.vert) ;

% 3. CALCULATION
%%%%%%%%%%%%%%%%

disp('flag 3.1.1')
tmp1 = C21/C11st ;
disp('flag 3.1.2')
tmp2 = C12/C22st ;
disp('flag 3.1.3')
tmp3 = C32/(- tmp1 * C12 + C22st ) ;
disp('flag 3.1.4')
tmp4 = C31/(- tmp2 * C21 + C11st ) ;
clear C21 C12 C11st C22st C32 C31

disp('flag 3.2.1')
tmp5 = tmp3*tmp1-tmp4 ;
tmp6 = tmp4*tmp2-tmp3 ;

disp('flag 3.2.2')
gama1 = - C33st - tmp5*C13 - tmp6*C23 ;
clear C33st C23 C13 tmp1 tmp2 tmp3 tmp4

disp('flag 3.3.1')
weight = 1/(2*pi*(sig1+sig2)) ;
G1 = bem_Gi_vert(tsph1.vert,dsph.dip,weight)' ;
weight = 1/(2*pi*(sig2+sig3)) ;
G2 = bem_Gi_vert(tsph2.vert,dsph.dip,weight)' ;
weight = 1/(2*pi*sig3) ;
G3 = bem_Gi_vert(tsph3.vert,dsph.dip,weight)' ;

disp('flag 3.3.2')
gama2 = tmp5*G1 + tmp6*G2 + G3 ;

%save self_infl tmp5 tmp6 gama1
clear G1 G2 G3 tmp5 tmp6

disp('flag 3.3.3')
L = gama1\gama2 ;

clear gama1 gama2
elapsed_time = etime(clock, TICTOC)

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

