function [tsurf_out] = split_tri(tsurf_in,option,Pvol)
%
% function [tsurf_out] = split_tri(tsurf_in,option) :
% refine mesh by dividing each triangle into 4 sub_triangle
% 	IN : tsurf_in.vert .tri .edges,option,Pvol
%	OUT : tsurf_out.split.vert .tri .neighb .nr
%
% structure of tsurf_out.split.neighb (TRANSPOSED !)
%
%    # of original pts -->
%  +---------------------------------------------+
%  | # of neighbours                             |
%  +---------------------------------------------+
%  |                                             |
%  |                                             |
%  | index of neighbours in tsurf_out.split.vert |
%  |   (by column)                               |
%  |                                             |
%  +---------------------------------------------+
%
% Split the triangles :
%
%	        1
%	        o
%	       / \		1,2,3 original vertices
%	      / I \		a,b,c new vertices
%	   a x-----x c		I,II,III,IV new triangles
%	    /\ IV  /\			I ->   1,a,c
%	   /  \   /  \			II ->  2,b,a
%	  / II \ / III\			III -> 3,c,b
%	 o------x------o		IV ->  a,b,c
%	2       b       3
%
% if 1 arg -> sphere 
% 	option=0
% if option=1
%	just take the midle of the edges
% elseif option=2
%	project midle edge on the true surface (scalp and brain)
% elseif option=3
%	project midle edge on the true surface for upper part only (skull),
% 	for lower part, take the middle of the edge
%


if nargin==0
	error('[tsurf_out] = split_tri(tsurf_in,option)') ;
elseif ~isfield(tsurf_in,'edges')
	error('I need the edges informations !  Use "measure_edges".') ;
end

vert = tsurf_in.vert ;
tri = tsurf_in.tri ;
Ntri = tsurf_in.nr(2) ;
Npt = tsurf_in.nr(1) ;
edge = tsurf_in.edges.ed ;
Ned = tsurf_in.edges.nr(1) ;

if nargin==1
	option=0 ;
	r_sph = norm(vert(1,:)) ;
elseif option==1
	disp('Middle of edges')
elseif (option==2 | option==3) 
	disp('Use true surface')
	if nargin==2
		Pvol = spm_get(1,'_bin.img','Bin image to use') ;
	end
	[DIM VOX SCALE TYPE OFFSET ORIGIN DESCRIP] = spm_hread(Pvol) ;
	Vv = spm_map(Pvol) ;
else
	error('Wrong option specification') ;
end

% parameters for the image Pvol
trsh_vol = .9 ; 	% for the projection of pts on the bin surf
n_div = 10 ;
w_inter = 15 ; 	% width of interval [-w w] in voxel
dr = 1/n_div ;
d_li = -w_inter:dr:w_inter ;
nd_li = length(d_li) ;
unit_li = ones(1,nd_li) ;
hold = 1 ;
load defin

% To be created !
vert_fine = zeros(Ned,3) ;
tri_fine = zeros(4*Ntri,3) ;
Max_neighb = 10 ;
neighbour = zeros(Npt,Max_neighb+1) ;

nrefin_pt = 1 ;

for i = 1:Ned
	ipt1 = edge(i,2) ;
	ipt2 = edge(i,3) ;
	new_pt = (vert(ipt1,:) + vert(ipt2,:))/2 ;

	if option==0
	% Sphere, take the middle and adjust for radius.
		new_pt = new_pt/norm(new_pt)*r_sph ;

	elseif option==1
	% Just keep the middle of segment
		% Nothing to do.
	elseif option==2
	% Project point on the surface
		itriI = edge(i,4) ;
		itriII = edge(i,5) ;
		iptI1 = tri(itriI,1) ;
		iptI2 = tri(itriI,2) ;
		iptI3 = tri(itriI,3) ;
		iptII1 = tri(itriII,1) ;
		iptII2 = tri(itriII,2) ;
		iptII3 = tri(itriII,3) ;
		ptI1 = vert(iptI1,:) ; 
		ptI2 = vert(iptI2,:) ; 
		ptI3 = vert(iptI3,:) ;
		ptII1 = vert(iptII1,:) ;
		ptII2 = vert(iptII2,:) ;
		ptII3 = vert(iptII3,:) ;

		c12 = ptI2 - ptI1 ;
		c13 = ptI3 - ptI1 ;
		normal_triI = [ c12(2)*c13(3)-c12(3)*c13(2) ; ...
				c12(3)*c13(1)-c12(1)*c13(3) ; ...
				c12(1)*c13(2)-c12(2)*c13(1) ] ;
		normal_triI = normal_triI/norm(normal_triI) ;
		c12 = ptII2 - ptII1 ;
		c13 = ptII3 - ptII1 ;
		normal_triII = [ c12(2)*c13(3)-c12(3)*c13(2) ; ...
				c12(3)*c13(1)-c12(1)*c13(3) ; ...
				c12(1)*c13(2)-c12(2)*c13(1) ] ;
		normal_triII = normal_triII/norm(normal_triII) ;

		normal_edge = (normal_triI+normal_triII) ...
				/ norm(normal_triI+normal_triII) ;

		line = new_pt'*unit_li + normal_edge*d_li ;
		val_line = spm_sample_vol(Vv,line(1,:),line(2,:),line(3,:), ...
						hold)/SCALE ;
%		pos = max(find(val_line>trsh_vol)) ;
		diff_val = diff(val_line) ;
		pos = min(find(diff_val<0)) +1 ;

		if ~(isempty(pos) | line(3,pos)<2)
			new_pt = line(:,pos)' ;
		end

	elseif option==3
	% Project point on the surface for upper part only,
	% for lower part, take the middle of the edge
		itriI = edge(i,4) ;
		itriII = edge(i,5) ;
		iptI1 = tri(itriI,1) ;
		iptI2 = tri(itriI,2) ;
		iptI3 = tri(itriI,3) ;
		iptII1 = tri(itriII,1) ;
		iptII2 = tri(itriII,2) ;
		iptII3 = tri(itriII,3) ;
		ptI1 = vert(iptI1,:) ; 
		ptI2 = vert(iptI2,:) ; 
		ptI3 = vert(iptI3,:) ;
		ptII1 = vert(iptII1,:) ;
		ptII2 = vert(iptII2,:) ;
		ptII3 = vert(iptII3,:) ;

		c12 = ptI2 - ptI1 ;
		c13 = ptI3 - ptI1 ;
		normal_triI = [ c12(2)*c13(3)-c12(3)*c13(2) ; ...
				c12(3)*c13(1)-c12(1)*c13(3) ; ...
				c12(1)*c13(2)-c12(2)*c13(1) ] ;
		normal_triI = normal_triI/norm(normal_triI) ;
		c12 = ptII2 - ptII1 ;
		c13 = ptII3 - ptII1 ;
		normal_triII = [ c12(2)*c13(3)-c12(3)*c13(2) ; ...
				c12(3)*c13(1)-c12(1)*c13(3) ; ...
				c12(1)*c13(2)-c12(2)*c13(1) ] ;
		normal_triII = normal_triII/norm(normal_triII) ;

		normal_edge = (normal_triI+normal_triII) ...
				/ norm(normal_triI+normal_triII) ;

		line = new_pt'*unit_li + normal_edge*d_li ;
		val_line = spm_sample_vol(Vv,line(1,:),line(2,:),line(3,:), ...
						hold)/SCALE;
%		pos = max(find(val_line>trsh_vol)) ;
		diff_val = diff(val_line) ;
		pos = min(find(diff_val<0)) ;
		if ((~isempty(pos)) & (line(3,pos+1)>(sl_skull+.5)))
			new_pt = line(:,pos+1)' ;
		end
	end
	vert_fine(i,:) = new_pt ;

	Nneigh_pt1 = neighbour(ipt1,1) ;
	Nneigh_pt2 = neighbour(ipt2,1) ;


	neighbour(ipt1,Nneigh_pt1+2) = i+Npt ;
	neighbour(ipt1,1) = Nneigh_pt1+1 ;
	neighbour(ipt2,Nneigh_pt2+2) = i+Npt ;
	neighbour(ipt2,1) = Nneigh_pt2+1 ;
end

Mneighb = max(neighbour(:,1)) ;
neighbour = neighbour(:,1:1+Mneighb) ;

% Organize triangles
for i = 1:Ntri
	ltri = find((edge(:,4)==i) | (edge(:,5)==i)) ;
	eda_tri = [tri(i,1) tri(i,2)] ;
	edb_tri = [tri(i,2) tri(i,3)] ;
	edc_tri = [tri(i,3) tri(i,1)] ;

	ipta = find(((eda_tri(1)==edge(ltri,2))&(eda_tri(2)==edge(ltri,3)))...
 		| ((eda_tri(1)==edge(ltri,3))&(eda_tri(2)==edge(ltri,2)))) ;
	iptb = find(((edb_tri(1)==edge(ltri,2))&(edb_tri(2)==edge(ltri,3)))...
 		| ((edb_tri(1)==edge(ltri,3))&(edb_tri(2)==edge(ltri,2)))) ;
	iptc = find(((edc_tri(1)==edge(ltri,2))&(edc_tri(2)==edge(ltri,3)))...
 		| ((edc_tri(1)==edge(ltri,3))&(edc_tri(2)==edge(ltri,2)))) ;
	pta = ltri(ipta) + Npt ;
	ptb = ltri(iptb) + Npt ;
	ptc = ltri(iptc) + Npt ;

	tri_fine((i-1)*4+1:i*4,:) = [tri(i,1) pta ptc ; ...
					 tri(i,2) ptb pta ; ...
					 tri(i,3) ptc ptb ; ...
					 pta ptb ptc] ;
end

tsurf_out = tsurf_in ;
tsurf_out.split.vert = [ vert ; vert_fine ] ;
tsurf_out.split.tri = tri_fine ;
tsurf_out.split.neighb = neighbour ;
tsurf_out.split.nr = [Npt+Ned 4*Ntri] ;
tsurf_out.info = [tsurf_in.info, ', split triangles'] ;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 


return


%        1
%        o
%       / \
%      / I \
%   a x-----x c
%    /\ IV  /\
%   /  \   /  \
%  / II \ / III\
% o------x------o
%2       b       3


