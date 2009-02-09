function Lel = elec_lin(electr,head_model,dipoles)

% Lel = bem_elec_lin(electr,[head_model,dipoles])
%
% Calculate the lead field for the specified electrodes
% uses the standard Cii, Cij, Gi function with the lin approximation.

TICTOC = clock ;

% 1. INITIALISATION
%%%%%%%%%%%%%%%%%%%
disp('flag 1')

load defin
if nargin==0
	error(['You must at least specify the electrodes !'])
elseif nargin==1
	n1 = n ; n2 = n ; n3 = n ;
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
elseif nargin==3
	% With real head model !
else
	error(['Wrong inputs : L = elec_lin2(electr,head_model,dipoles) ;'])
end

Ntot_vert = tsph1.nr(1) + tsph2.nr(1) + tsph3.nr(1)
defl = 1/Ntot_vert

Nel  = electr.nr(1) ;
El_vert = tsph3.tri(electr.ctri,:) ;
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

% 2. MAIN MATRICES, P1
%%%%%%%%%%%%%%%%%%%%%%

disp('flag 2.2')
weight = (sig1-sig2)/((sig1+sig2)*2*pi) ;
C11st = bem_Cii_lin(tsph1.tri,tsph1.vert,weight,defl,tsph1.pt4.vert) ;
weight = (sig2-sig3)/((sig1+sig2)*2*pi) ;
C12 = bem_Cij_lin(tsph1.vert,tsph2.vert,tsph2.tri,weight,defl) ;

disp('flag 2.3')
weight = (sig1-sig2)/((sig2+sig3)*2*pi) ;
C21 = bem_Cij_lin(tsph2.vert,tsph1.vert,tsph1.tri,weight,defl) ;
weight = (sig2-sig3)/((sig2+sig3)*2*pi) ;
C22st = bem_Cii_lin(tsph2.tri,tsph2.vert,weight,defl,tsph2.pt4.vert) ;

disp('flag 2.4')
weight = (sig1-sig2)/(sig3*2*pi) ;
C31 = Perm_r*bem_Cij_lin(tsph3.vert,tsph1.vert,tsph1.tri,weight,defl) ;
weight = (sig2-sig3)/(sig3*2*pi) ;
C32 = Perm_r*bem_Cij_lin(tsph3.vert,tsph2.vert,tsph2.tri,weight,defl) ;

% 3. CALCULATION, P1
%%%%%%%%%%%%%%%%%%%%

disp('flag 3.1.1')
tmp1 = C21/C11st ;
tmp2 = C12/C22st ;
disp('flag 3.1.2')
tmp3 = C32/(- tmp1 * C12 + C22st ) ;
tmp4 = C31/(- tmp2 * C21 + C11st ) ;
clear C21 C12 C11st C22st C32 C31

disp('flag 3.1.3')
tmp5 = tmp3*tmp1-tmp4 ;
tmp6 = tmp4*tmp2-tmp3 ;

% 4. MAIN MATRICES, P2
%%%%%%%%%%%%%%%%%%%%%%

disp('flag 4')
weight = sig3/((sig1+sig2)*2*pi) ;
C13 = bem_Cij_lin(tsph1.vert,tsph3.vert,tsph3.tri,weight,defl)*Perm_c ;
weight = sig3/((sig2+sig3)*2*pi) ;
C23 = bem_Cij_lin(tsph2.vert,tsph3.vert,tsph3.tri,weight,defl)*Perm_c ;
weight = 1/(2*pi) ;
C33st = Perm_r*bem_Cii_lin(tsph3.tri,tsph3.vert,weight,defl,tsph3.pt4.vert)*Perm_c ;


% 5. CALCULATION, P2
%%%%%%%%%%%%%%%%%%%%

disp('flag 5.1')
gama1 = - C33st - tmp5*C13 - tmp6*C23 ;
clear C33st C23 C13 tmp1 tmp2 tmp3 tmp4

A = gama1(1:nRem_vert,1:nRem_vert) ;
B = gama1(1:nRem_vert,nRem_vert+1:tsph3.nr(1)) ;
C = gama1(nRem_vert+1:tsph3.nr(1),1:nRem_vert) ;
D = gama1(nRem_vert+1:tsph3.nr(1),nRem_vert+1:tsph3.nr(1)) ;
clear gama1

disp('flag 5.2')
CAi = C/A ;
clear C A

F = D - CAi*B ;
Fi = inv(F) ;
gama_el = [-Fi*CAi Fi] ;

disp('flag 5.3')
IFS1 = gama_el*tmp5 ; % Intermediate Forward Solution
IFS2 = gama_el*tmp6 ;
IFS3 = gama_el ;
clear D B tmp5 tmp6 gama_el

save IFS IFS1 IFS2 IFS3

% 6. MAIN MATRICES, P3
%%%%%%%%%%%%%%%%%%%%%%

disp('flag 6.1')
weight = 1/(2*pi*(sig1+sig2)) ;
G1 = bem_Gi_vert(tsph1.vert,dsph.dip,weight)' ;
weight = 1/(2*pi*(sig2+sig3)) ;
G2 = bem_Gi_vert(tsph2.vert,dsph.dip,weight)' ;
weight = 1/(2*pi*sig3) ;
G3 = Perm_r * bem_Gi_vert(tsph3.vert,dsph.dip,weight)' ;

Lvert = (IFS1*G1) + (IFS2*G2) + (IFS3*G3) ;

Lel =  zeros(Nel,3*dsph.nr(1)) ;
for i=1:Nel
	r1 = tsph3.vert(tsph3.tri(electr.ctri(i),1),:) ;
	r2 = tsph3.vert(tsph3.tri(electr.ctri(i),2),:) ;
	r3 = tsph3.vert(tsph3.tri(electr.ctri(i),3),:) ;
	z1 = [r2(2)*r3(3)-r3(2)*r2(3) ; r3(1)*r2(3)-r2(1)*r3(3) ; ...
		r2(1)*r3(2)-r3(1)*r2(2)] ; 
	z2 = [r3(2)*r1(3)-r1(2)*r3(3) ; r1(1)*r3(3)-r3(1)*r1(3) ; ...
		r3(1)*r1(2)-r1(1)*r3(2)] ;
	z3 = [r1(2)*r2(3)-r2(2)*r1(3) ; r2(1)*r1(3)-r1(1)*r2(3) ; ...
		r1(1)*r2(2)-r2(1)*r1(2)] ;
	h1 = electr.coord(i,:)*z1/(r1*z1) ;
	h2 = electr.coord(i,:)*z2/(r2*z2) ;
	h3 = electr.coord(i,:)*z3/(r3*z3) ;
	Lel(i,:) = Lvert(3*i-2,:)*h1 + Lvert(3*i-1,:)*h2 + Lvert(3*i,:)*h3 ;
end

save Lel Lel

elapsed_time = etime(clock, TICTOC)

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
