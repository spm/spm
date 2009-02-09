function tsurf_out = tri_4thpt(tsurf_in,option,Pvol)
% function tsurf_out = tri_4thpt(tsurf_in,option) :
%
% for each triangle, calculate the intersection of the mediatrice and 
% project it on the real surface.
%	IN : tsurf_in.vert .tri, option, Pvol
%	OUT : tsurf_out.pt4.vert
%
%	tsurf_out.pt4 (Ntri x 3) : coordinates of the center of tri.
%
% if 1 arg or option=0 -> sphere 
% 	
% elseif option=1 -> brain or scalp true surface 
%
% elseif option=2 -> skull : up, true surface ; low, local sphere

if nargin==1
	option = 0 ;
end

Ntri = tsurf_in.nr(2) ;
tri = tsurf_in.tri ;
vert = tsurf_in.vert' ;
load defin

if option==0
	r_sph = norm(vert(:,1)) ;
elseif ((option==1) | (option==2))
	% brain, skull or scalp bin file needed.
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
w_inter = 14 ; 	% width of interval [-w w] in voxel
dr = 1/n_div ;
d_li = -w_inter:dr:w_inter ;
nd_li = length(d_li) ;
unit_li = ones(1,nd_li) ;
hold = 1 ;


tr_4pt = zeros(Ntri,3) ;
for i = 1:Ntri
	ipt1 = tri(i,1) ;
	ipt2 = tri(i,2) ;
	ipt3 = tri(i,3) ;
	pt1 = vert(:,ipt1) ; % vertical vector,
	pt2 = vert(:,ipt2) ; % because vert has been transposed
	pt3 = vert(:,ipt3) ;

	c12 = pt2 - pt1 ;
	c13 = pt3 - pt1 ;
	normal_tri = [ c12(2)*c13(3)-c12(3)*c13(2) ; ...
			c12(3)*c13(1)-c12(1)*c13(3) ; ...
			c12(1)*c13(2)-c12(2)*c13(1) ] ;
	normal_tri = normal_tri/norm(normal_tri) ;

	pt_c = (pt1+pt2+pt3)/3 ;

% Projection on the surface
	if option==0
	% Sphere, adjust for radius.
		tr_4pt(i,:) = pt_c'/norm(pt_c)*r_sph ;

	elseif option==1
	% brain and scalp true surface, except for z<2 (scalp bootom)
		line = pt_c*unit_li + normal_tri*d_li ;
		val_line = spm_sample_vol(Vv,line(1,:),line(2,:),line(3,:), ...
						hold)/SCALE;
%		pos = max(find(val_line>trsh_vol)) ;
		diff_val = diff(val_line) ;
		pos = min(find(diff_val<0)) +1 ;

		if (isempty(pos) | line(3,pos)<2)
			tr_4pt(i,:) = approx_sph(i,ipt1,ipt2,ipt3, ...
					pt1,pt2,pt3,pt_c,tri,vert) ;
		else
			tr_4pt(i,:) = line(:,pos)' ;
		end


	elseif option==2
	% Skull : use true surface for upper part, local sphere for lower part.
		line = pt_c*unit_li + normal_tri*d_li ;
		val_line = spm_sample_vol(Vv,line(1,:),line(2,:),line(3,:), ...
						hold)/SCALE;
%		pos = max(find(val_line>trsh_vol)) ;
		diff_val = diff(val_line) ;
		pos = min(find(diff_val<0)) ;

		if ((~isempty(pos)) & (line(3,pos+1)>(sl_skull+.5)))
		% use true surface
			tr_4pt(i,:) = line(:,pos+1)' ;
		else
		% find local sphere and project pt_c on it.
			tr_4pt(i,:) = approx_sph(i,ipt1,ipt2,ipt3, ...
					pt1,pt2,pt3,pt_c,tri,vert) ;
		end
	end
end

tsurf_out = tsurf_in ;
tsurf_out.pt4.vert = tr_4pt ;
tsurf_out.info = [tsurf_in.info, ', triangle 4th pt'] ;

if ((option==1) | (option==2))
	spm_unmap(Vv)
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sph_4thpt = approx_sph(i,ipt1,ipt2,ipt3,pt1,pt2,pt3,pt_c,tri,vert)

% Need a 4th point to define the sphere : take the closest pt to pt1,2 or 3.
% For pt1 : find tri containing ipt1, discard tri i
%	    -> list of close points, discard ipt1,2,3 and doubles
%	    -> measure distance between close points and pt1
% do the same for pt2 and pt3
% Compare the distances, take the smallest -> 4th point use to find local sph.

closest=[] ; dist_closest = [] ;
for j=1:3
% run through the 3 points ipt1,2,3
	eval(['ipt = ipt',num2str(j),';']) ;
	eval(['pt = pt',num2str(j),';']) ;
	li_tri = find( (tri(:,1)==ipt) | (tri(:,2)==ipt) | (tri(:,3)==ipt) ) ;
	li_tri(find(li_tri==i))=[] ;
	li_pts = tri(li_tri,:) ; li_pts = li_pts(:) ;
	li_pts(find((li_pts==ipt1) | (li_pts==ipt2) | (li_pts==ipt3)))=[] ;
	for k=1:length(closest)
		li_pts(find(li_pts==closest(k))) = [] ;
	end
	% Remove doubles
	if ~isempty(li_pts)
		li_pts = sort(li_pts) ;
		d = diff(li_pts) ;
		li_pts = li_pts([find(d~=0) ; length(li_pts)]) ;
		Npts = length(li_pts) ;
		% Measure dist and keep the min
		dist = sqrt(sum((vert(:,li_pts) - pt*ones(1,Npts)).^2)) ;
%		mdist = find(dist==min(dist)) ;
%		closest = [closest li_pts(mdist)'] ;
%		dist_closest = [dist_closest dist(mdist)] ;
		closest = [closest li_pts'] ;
		dist_closest = [dist_closest dist] ;
	end
end
[dist_closest,perm] = sort(dist_closest) ;
closest = closest(perm) ;

pt4 = vert(:,closest(1)) ;
rA = 2 ;
while rA==2
% to avoid pb of 4 points on 1 plane, rA should be 3
	% Generate the sphere going through the vertices 
	% of the triangles + the closest point.
	% cfr paper by L. Heller
	m2 = (pt1+pt2)/2 ; m3 = (pt1+pt3)/2 ; m4 = (pt1+pt4)/2 ; 
	d2 = pt1-pt2 ; d3 = pt1-pt3 ; d4 = pt1-pt4 ;
	A = [ d2' ; d3' ; d4' ] ;
	b = [ m2'*d2 ; m3'*d3 ; m4'*d4 ] ;
	rA = rank(A)	 ;
	if rA~=3	
		closest	(1)=[] ;
		dist_closest(1)=[] ;
		if isempty(closest)
			c12 = pt2 - pt1 ;
			c13 = pt3 - pt1 ;
			normal_tri = [ c12(2)*c13(3)-c12(3)*c13(2) ; ...
				c12(3)*c13(1)-c12(1)*c13(3) ; ...
				c12(1)*c13(2)-c12(2)*c13(1) ] ;
			normal_tri = normal_tri/norm(normal_tri) ;
			pt4 = pt_c + normal_tri*.001 ;
			% No other solution than taking the pt4 
			% a bit above the triangle...
		else
			pt4 = vert(:,closest(1)) ;
		end
	end
end
centre_sph = A\b ;
R = norm(pt1-centre_sph) ;
ve = pt_c-centre_sph ;
ve = ve/norm(ve) ;

sph_4thpt = (centre_sph + R*ve)' ;

