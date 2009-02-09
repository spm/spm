function tsurf_out = measure_edges(tsurf_in)
% 
% function tsurf_out = measure_edges(tsurf_in)
%
% 	IN : tsurf_in.vert .tri .nr
%	OU : tsurf_out.edges.ed .nr
%
% 	tsurf_out.edges.ed :
% 		ed(:,1) = length,
%		ed(:,[2 3]) = indexes of 2 points,
%		ed(:,[4 5]) = indexes of 2 triangles.
%	tsurf_out.edges.nr :
%		nr(1) : number of edges,
%		nr([2 3]) : mean and std of edges length.

if nargin==0
	error(['tsurf_out = measure_edges(tsurf_in)']) ;
end

Ntri = tsurf_in.nr(2) ;
ind_tr = tsurf_in.tri ;
xyz = tsurf_in.vert ;

Ned = 3*Ntri/2 ;
edges = zeros(Ned,5) ;
i_ed = 1 ;

% Go through all the triangles
for i=1:Ntri
%	if rem(i,400)==0 i , end ;
	% 3 edges of the ith triangle (counter-clock wise)
	e1 = [ind_tr(i,1) ind_tr(i,2)] ; 
	e2 = [ind_tr(i,2) ind_tr(i,3)] ;
	e3 = [ind_tr(i,3) ind_tr(i,1)] ;
	% check if edge was already seen.
	p1 = find((edges(:,2)==e1(2)) & (edges(:,3)==e1(1))) ;
	p2 = find((edges(:,2)==e2(2)) & (edges(:,3)==e2(1))) ;
	p3 = find((edges(:,2)==e3(2)) & (edges(:,3)==e3(1))) ;
	if isempty(p1)
		% new edge
		le1 = norm(xyz(e1(1),:)-xyz(e1(2),:)) ;
		edges(i_ed,:) = [le1 e1 i 0] ;
		i_ed = i_ed+1 ;
	else
		% 2nd triangle for this edge
		edges(p1,5) = i ;
	end
	if isempty(p2)
		le2 = norm(xyz(e2(1),:)-xyz(e2(2),:)) ;
		edges(i_ed,:) = [le2 e2 i 0] ;
		i_ed = i_ed+1 ;
	else
		edges(p2,5) = i ;
	end
	if isempty(p3)
		le3 = norm(xyz(e3(1),:)-xyz(e3(2),:)) ;
		edges(i_ed,:) = [le3 e3 i 0] ;
		i_ed = i_ed+1 ;
	else
		edges(p3,5) = i ;
	end
end
me_len  = mean(edges(:,1)) ;
std_len = std(edges(:,1)) ;

tsurf_out = tsurf_in ;
tsurf_out.edges.ed = edges ;
tsurf_out.edges.nr = [Ned me_len std_len] ;
tsurf_out.info = [tsurf_in.info, ', edge measured'] ;

