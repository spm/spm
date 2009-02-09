%%
% 1. Head model creation: meshed surfaces, sources and electrodes
%%%%%%%%%%%%%%%%%%%%%%%%%

% 1.1 Create surfaces
%====================
% create unit sphere
[pnt, tri] = icosahedron162;

vol = [];
vol.cond = [1 1/80 1];
vol.source = 1; % index of source compartment
vol.skin   = 3; % index of skin surface
% brain
vol.bnd(1).pnt = pnt*88;
vol.bnd(1).tri = tri;
% skull
vol.bnd(2).pnt = pnt*92;
vol.bnd(2).tri = tri;
% skin
vol.bnd(3).pnt = pnt*100;
vol.bnd(3).tri = tri;

% Build Triangle 4th point
vol = triangle4pt(vol); % TO BE WRITTEN
% vol.bnd(1-2-3).pnt4, i.e. One pt4 per triangle !

% 1.2 Create electrodes
%======================
pnt = randn(64,3);
for i=1:64
  elec.pnt(i,:) = pnt(i,:) / norm(pnt(i,:)); % scale towards the skin surface
  elec.label{i} = sprintf('%02d', i);
end

% 1.3 create sources
%===================

Rbr = norm(vol.bnd(1).pnt(1,:));
tmp = (rand(4000,3)-.5)*Rbr;
dip.pnt = tmp(sqrt(sum(tmp.^2,2))<Rbr*.95,:);
    % Keep those inside the brain volume

%%
% 2. BEM model estimation, only for the scalp surface
%========================
% see below for all the surfaces.

defl =[ 0 0 1/size(vol.bnd(vol.skin).pnt,1)] ;
% ensure deflation for skin surface, i.e. average reference over skin

% NOTE:
% Calculation proceeds by estimating each submatrix C_ij and combine them. 
% There are 2 options:
% - calculating the matrices once, as it takes some time, keep them in
%   memory and use them the 2-3 times they're needed.
% - calculating the matrices every time they're needed, i.e. 2-3 times
% The former option is faster but requires more memory space as up to *8* 
% square matrices of size C_ij have to be kept in memory at once.
% The latter option requires less memory, but would take much more time to
% estimate.
% This faster but memory hungry solution is implemented here.

% Deal first with surface 1 and 2 (inner and outer skull
%--------------------------------

% NOTE:
% C11st/C22st/C33st are simply the matrix C11/C22/C33 minus the identity
% matrix, i.e. C11st = C11-eye(N)

weight = (vol.cond(1)-vol.cond(2))/((vol.cond(1)+vol.cond(2))*2*pi);
C11st = bem_Cii_lin(vol.bnd(1).tri,vol.bnd(1).pnt, ...
                weight,defl(1),vol.bnd(1).pnt4);
weight = (vol.cond(1)-vol.cond(2))/((vol.cond(2)+vol.cond(3))*2*pi) ;
C21 = bem_Cij_lin(vol.bnd(2).pnt,vol.bnd(1).pnt,vol.bnd(1).tri, ...
              weight,defl(1)) ;
tmp1 = C21/C11st ;

weight = (vol.cond(2)-vol.cond(3))/((vol.cond(1)+vol.cond(2))*2*pi);
C12 = bem_Cij_lin(vol.bnd(1).pnt,vol.bnd(2).pnt,vol.bnd(2).tri, ...
              weight,defl(2)) ;
weight = (vol.cond(2)-vol.cond(3))/((vol.cond(2)+vol.cond(3))*2*pi) ;
C22st = bem_Cii_lin(vol.bnd(2).tri,vol.bnd(2).pnt, ...
                weight,defl(2),vol.bnd(2).pnt4) ;
tmp2 = C12/C22st ;

% Combine with the effect of surface 3 (scalp) on the first 2
%------------------------------------------------------------
weight = (vol.cond(1)-vol.cond(2))/(vol.cond(3)*2*pi) ;
C31 = bem_Cij_lin(vol.bnd(3).pnt,vol.bnd(1).pnt,vol.bnd(1).tri, ...
                weight,defl(1)) ;
tmp4 = C31/(- tmp2 * C21 + C11st ) ;
clear C31 C21 C11st
            
weight = (vol.cond(2)-vol.cond(3))/(vol.cond(3)*2*pi) ;
C32 = bem_Cij_lin(vol.bnd(3).pnt,vol.bnd(2).pnt,vol.bnd(2).tri, ...
                weight,defl(2)) ;
tmp3 = C32/(- tmp1 * C12 + C22st ) ;
clear  C12 C22st C32

tmp5 = tmp3*tmp1-tmp4 ;
tmp6 = tmp4*tmp2-tmp3 ;
clear tmp1 tmp2 tmp3 tmp4

% Finally include effect of surface 3 on the others
%--------------------------------------------------
% As the gama1 intermediate matrix is built as the sum of 3 matrices, I can
% spare some memory by building them one at a time, and summing directly
weight = vol.cond(3)/((vol.cond(1)+vol.cond(2))*2*pi) ;
Ci3 = bem_Cij_lin(vol.bnd(1).pnt,vol.bnd(3).pnt,vol.bnd(3).tri, ...
                weight,defl(3)) ;
gama1 = - tmp5*Ci3 ; % gama1 = - tmp5*C13;

weight = vol.cond(3)/((vol.cond(2)+vol.cond(3))*2*pi) ;
Ci3 = bem_Cij_lin(vol.bnd(2).pnt,vol.bnd(3).pnt,vol.bnd(3).tri, ...
                weight,defl(3)) ;
gama1 = gama1 - tmp6*Ci3; % gama1 = - tmp5*C13 - tmp6*C23;

weight = 1/(2*pi) ;
Ci3 = bem_Cii_lin(vol.bnd(3).tri,vol.bnd(3).pnt, ...
                weight,defl(3),vol.bnd(3).pnt4) ;
gama1 = gama1 - Ci3; % gama1 = - tmp5*C13 - tmp6*C23 - C33st ;
clear Ci3 tmp1 tmp2 tmp3 tmp4

% Other wise I can create all 3 C_i3 matrices at once and combine them  
% all together at the very end, as here under
weight = vol.cond(3)/((vol.cond(1)+vol.cond(2))*2*pi) ;
C13 = bem_Cij_lin(vol.bnd(1).pnt,vol.bnd(3).pnt,vol.bnd(3).tri, ...
                weight,defl(3)) ;
weight = vol.cond(3)/((vol.cond(2)+vol.cond(3))*2*pi) ;
C23 = bem_Cij_lin(vol.bnd(2).pnt,vol.bnd(3).pnt,vol.bnd(3).tri, ...
                weight,defl(3)) ;
weight = 1/(2*pi) ;
C33st = bem_Cii_lin(vol.bnd(3).tri,vol.bnd(3).pnt, ...
                weight,defl(3),vol.bnd(3).pnt4) ;

gama1 = - tmp5*C13 - tmp6*C23 - C33st ;
clear C33st C23 C13 tmp1 tmp2 tmp3 tmp4
% 
%%
% % 3. Build the leadfield from the BEM and dipoles, NOT as done in FT!!!
% %================================================
% % For FT implementation, use code below
% 
% % Direct leadfield matrices, for the 3 surfaces.
% %--------------------------
% % Build the direct leadfield for each model surface, once at a time and
% % combine it on the fly, saving some memory.
% 
% % weight = 1/(2*pi*(vol.cond(1)+vol.cond(2))) ;
% % G = bem_Gi_vert(vol.bnd(1).pnt,dip.pnt,weight)' ;
% co = (vol.cond(1)+vol.cond(2))/2 ;
% G = inf_medium_leadfield(dip.pnt, vol.bnd(1).pnt, co) ;
% gama2 = tmp5*G;
% % clear G tmp5
% 
% % weight = 1/(2*pi*(vol.cond(2)+vol.cond(3))) ;
% % G = bem_Gi_vert(vol.bnd(2).pnt,dip.pnt,weight)' ;
% co = (vol.cond(2)+vol.cond(3))/2 ;
% G = inf_medium_leadfield(dip.pnt, vol.bnd(2).pnt, co) ;
% gama2 = gama2 + tmp6*G;
% % clear G tmp6
% 
% % weight = 1/(2*pi*vol.cond(3)) ;
% % G = bem_Gi_vert(vol.bnd(3).pnt,dip.pnt,weight)' ;
% co = vol.cond(3)/2 ;
% G = inf_medium_leadfield(dip.pnt, vol.bnd(3).pnt, co) ;
% gama2 = gama2 + G ;
% % clear G 
% 
% % % Or build all 3 direct lf matrices and combine them all together at the
% % % end
% % weight = 1/(2*pi*(vol.cond(1)+vol.cond(2))) ;
% % G1 = bem_Gi_vert(vol.bnd(1).pnt,dip.pnt,weight)' ;
% % weight = 1/(2*pi*(vol.cond(2)+vol.cond(3))) ;
% % G2 = bem_Gi_vert(vol.bnd(2).pnt,dip.pnt,weight)' ;
% % weight = 1/(2*pi*vol.cond(3)) ;
% % G3 = bem_Gi_vert(vol.bnd(3).pnt,dip.pnt,weight)' ;
% % gama2 = tmp5*G1 + tmp6*G2 + G3 ;
% % clear G1 G2 G3 tmp5 tmp6
% 
% % Leadfiled for the scalp surface mesh !
% %---------------------------------------
% lf = gama1\gama2 ;
% 
% lf_3_3 = lf;

%%
% 4. Build leadfield and stuff as in FT !
%========================================

% Build system matrix
%--------------------
i_gama1 = inv(gama1);
vol.mat = [i_gama1*tmp5 i_gama1*tmp6 i_gama1];

% infinite medium leadfield
%--------------------------
cond = [vol.cond 0]; % add up the conductivity of air for simplicity
lf = [];
for ii=1:3
%     weight = 1/(2*pi*(cond(ii)+cond(ii+1))) ;
%     lf = [lf ; bem_Gi_vert(vol.bnd(ii).pnt,dip.pnt,weight)'] ;
    co = (cond(ii)+cond(ii+1))/2 ;
    lf = [lf ; inf_medium_leadfield(dip.pnt, vol.bnd(ii).pnt, co)] ;
end

% leadfield for bounded medium on scalp vertices
%-----------------------------------------------
lf = vol.mat * lf;
lf_4_3 = lf;

%%
% 5. BEM model estimation for all the vertices of the 3 surfaces
%    and leadfield estimation
%===============================================================

% deflation coeficient, ensuring zero mean scalp potential
defl =[ 0 0 1/size(vol.bnd(vol.skin).pnt,1)] ;
cond = [vol.cond 0]; % add up the conductivity of air for simplicity

% Note Cst = C-I, therefore system matrix = -Cst^-1
Cst = [] ;
tmp = cell(1,3) ;
for ii=1:3 % 3 lines
    for jj=1:3 % 3 columns
        if ii==jj
            weight = (cond(jj)-cond(jj+1)) / ((cond(jj)+cond(jj+1))*2*pi);
            tmp{jj} = bem_Cii_lin(vol.bnd(jj).tri,vol.bnd(jj).pnt, ...
                                weight,defl(jj),vol.bnd(jj).pnt4);
        else
            weight = (cond(jj)-cond(jj+1)) / ((cond(ii)+cond(ii+1))*2*pi);
            tmp{jj} = bem_Cij_lin(vol.bnd(ii).pnt,vol.bnd(jj).pnt, ...
                        vol.bnd(jj).tri, weight,defl(2)) ;
        end
    end
    Cst = [Cst ; tmp{1} tmp{2} tmp{3}];
end
vol.mat = -inv(Cst);

% infinite medium leadfield
%--------------------------
lf = [];
for ii=1:3
%     weight = 1/(2*pi*(cond(ii)+cond(ii+1))) ;
%     lf = [lf ; bem_Gi_vert(vol.bnd(ii).pnt,dip.pnt,weight)'] ;
    co = (cond(ii)+cond(ii+1))/2 ;
    lf = [lf ; inf_medium_leadfield(dip.pnt, vol.bnd(ii).pnt, co)] ;
end

% leadfield for bounded medium on all vertices
%---------------------------------------------
lf = vol.mat * lf;

lf_5 = lf; lf_5_3 = lf_5((162*2+1):end,:);


return

%%
% Check results
max(abs(lf_5_3(:)))
max(abs(lf_3_3(:)))
max(abs(lf_4_3(:)))

d_53 = lf_5_3-lf_3_3;
d_54 = lf_5_3-lf_4_3;
d_34 = lf_3_3-lf_4_3;

imat({d_53,d_54,d_34},2)

max(abs(d_53(:)))
max(abs(d_54(:)))
max(abs(d_34(:)))

% compare  my routine and Robert's for direct lf
weight = 1/(2*pi*(vol.cond(1)+vol.cond(2))) ;
tic
G1 = bem_Gi_vert(vol.bnd(1).pnt,dip.pnt,weight)' ;
toc

co = (vol.cond(1)+vol.cond(2))/2
tic
[lf] = inf_medium_leadfield(dip.pnt, vol.bnd(1).pnt, co);
toc

imat(G1-lf,2)
