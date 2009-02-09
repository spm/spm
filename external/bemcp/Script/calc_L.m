%function calc_L

% function calc_L
%
% Calculates the final solution : Lvert and Lel
% using the results of calc_IFS
%
% Loads IFS and param
% Saves Lvert and Lel

%load Param
%load Param_n40_d5271
load IFS

Lvert 	= zeros(nEl_vert,3*dsph.nr(1)) ;
Lel	= zeros(Nel,3*dsph.nr(1)) ;

we = weight([10 11 12]) ;

for i=1:n_shareG
	if i~=n_shareG
		disp([num2str(n_shareG-i),' parts left to process.']) ;
		eval(['dipole_loc = dsph.dip(share', num2str(i) ,',:) ;']) ;
		eval(['Lvert(:,(i-1)*shareG*3+1:i*shareG*3) = ',...
			'end_IFS(dipole_loc,head,IFS1,IFS2,IFS3,Perm_r,we) ;']);
	else
		disp('Processing last part.') ;
		eval(['dipole_loc = dsph.dip(share', num2str(i) ,',:) ;']) ;
		eval(['Lvert(:,(i-1)*shareG*3+1:dsph.nr(1)*3) = ',...
			'end_IFS(dipole_loc,head,IFS1,IFS2,IFS3,Perm_r,we) ;']);
	end
end

% save Lvert_ALL Lvert

for i=1:Nel
	r1 = tsph3.vert(El_vert(3*i-2),:) ;
	r2 = tsph3.vert(El_vert(3*i-1),:) ;
	r3 = tsph3.vert(El_vert(3*i),:) ;
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

save Lel_ALL5 Lel
