function tsph = spm_eeg_inv_TesSph(r,n);

%=======================================================================
% FORMAT tsph = spm_eeg_inv)TesSph(r,n);
%
% Generate a structure 'tsph' containing a tesselated sphere.
%
% Input : 
% r             - radius of the sphere
% n             - number of 'latitude' divisions on the sphere. It MUST be even!
%               ( Nvert = 5/4 * n^2 + 2 )
%
% Output : 
% tsph .vert    - vertices coordinates (3 x Nvert)
%      .tri     - triangle patches (3 x Ntri)
%	   .info    - info string
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips & Jeremie Mattout
% $Id: spm_eeg_inv_TesSph.m 716 2007-01-16 21:13:50Z karl $

if nargin == 0 
    n = 40;
	r = 1;
elseif nargin == 1
    n = 40;
elseif nargin == 2
    n = round( sqrt(4*(n - 2)/5) );
else
    error('Error: wrong input format')
end
	
[X_sp,Y_sp,Z_sp,npt_tr,npt_qr,npt_sph] = point(n,r) ;
[ind_tr] = tri(npt_tr,npt_qr,npt_sph) ;
	
tsph.vert = [X_sp ; Y_sp ; Z_sp] ;
tsph.tri = ind_tr ;
tsph.nr = [length(X_sp) size(ind_tr,2)] ;
tsph.info = ['Tes sph : n=',num2str(n),', r=',num2str(r)] ;
varargout{1} = tsph;


%________________________________________________________________________
% FORMAT [X_sp,Y_sp,Z_sp,npt_tr,npt_qr,npt_sph] = point(n,r)
% Generate a set of points on a sphere.
% n = the number of latitudes on the sphere, it MUST be even!
% r = radius of the sphere, 1 by default.
%------------------------------------------------------------------------
function [X_sp,Y_sp,Z_sp,npt_tr,npt_qr,npt_sph] = point(n,r)
d_r = pi/180; dth = 180/n*d_r ;
X_qr=[] ; Y_qr=[] ; Z_qr=[] ; X_sp=[] ; Y_sp=[] ; Z_sp=[] ;

% Coord of pts on a 5th of sphere WITHOUT the top and bottom points
for i=1:n/2         % Upper part
	the = i*dth ;
	cth = cos(the) ; sth = sin(the) ;
	dphi = 72/i*d_r ;
	for j=1:i
		phi = (j-1)*dphi ;
		cph = cos(phi) ; sph = sin(phi) ;
		X_qr = [X_qr sth*cph] ;
		Y_qr = [Y_qr sth*sph] ;
		Z_qr = [Z_qr cth] ;
	end
end
for i=(n/2+1):(n-1) % Lower part
	the = i*dth ;
	cth = cos(the) ; sth = sin(the) ;
	ii = n-i ;
	dphi = 72/ii*d_r ;
	for j=1:ii
		phi = (j-1)*dphi ;
		cph = cos(phi) ; sph = sin(phi) ;
		X_qr = [X_qr sth*cph] ;
		Y_qr = [Y_qr sth*sph] ;
		Z_qr = [Z_qr cth] ;
	end
end
% Each part is copied with a 72deg shift
dphi=-72*d_r ;
for i=0:4
	phi=i*dphi ;
	cph=cos(phi) ;
	sph=sin(phi) ;
	X_sp=[X_sp cph*X_qr+sph*Y_qr] ;
	Y_sp=[Y_sp -sph*X_qr+cph*Y_qr] ;
	Z_sp=[Z_sp Z_qr] ;
end
% add top and bottom points
X_sp=[0 X_sp 0]*r ;
Y_sp=[0 Y_sp 0]*r ;
Z_sp=[1 Z_sp -1]*r ;

npt_tr = [[0:n/2] [n/2-1:-1:0]]; % Nbr of points on each slice for a 5th of sph.
npt_qr = sum(npt_tr) ; % total nbr of pts per 5th of sphere
npt_sph=5/4*n^2+2 ; % or 5*npt_qr+2 ; % total nbr of pts on sph
return

%________________________________________________________________________
% FORMAT [ind_tr] = tri(npt_tr,npt_qr,npt_sph)
% Order the pts created by the function 'point' into triangles
% Each triangle is represented by 3 indices correponding to the XYZ_sp
%
% This works only because I know in wich order the vertices are generated
% in the 'point' function.
%------------------------------------------------------------------------
function [ind_tr] = tri(npt_tr,npt_qr,npt_sph)
n = sqrt(4/5*(npt_sph-2)) ;

% First 5th of sphere only
for i=0:(n/2-1)         % upper half
	if i==0		    	% 1st triangle at the top
		ind_qr=[1 2 2+npt_qr]' ;
	else		    	% other triangles
		for j=0:npt_tr(i+1)
			if j==npt_tr(i+1)-1			% last but 1 pt on slice
				x1=sum(npt_tr(1:i))+j+2 ;
				x12=x1+npt_tr(i+1) ;
				x13=x12+1 ;
				x22=x13 ;
				x23=sum(npt_tr(1:i))+2+npt_qr ;
				ind_qr=[ind_qr [x1 x12 x13]' [x1 x22 x23]'] ;
			elseif j==npt_tr(i+1)		% last pt on slice
				x1=sum(npt_tr(1:i))+2+npt_qr ;
				x2=sum(npt_tr(1:(i+1)))+j+2 ;
				x3=x1+npt_tr(i+1) ;
				ind_qr=[ind_qr [x1 x2 x3]'] ;
			else        			% other pts on slice
				x1=sum(npt_tr(1:i))+j+2 ;
				x12=x1+npt_tr(i+1) ;
				x13=x12+1 ;
				x22=x13 ;
				x23=x1+1 ;
				ind_qr=[ind_qr [x1 x12 x13]' [x1 x22 x23]'] ;
			end
		end
	end
end
for i=(n/2+1):n         % lower half
	if i==n	        	% last triangle at the bottom
		ind_qr=[ind_qr [5*npt_qr+2 2*npt_qr+1 npt_qr+1]' ] ;
	else            	% other triangles
		for j=0:npt_tr(i+1)
			if j==npt_tr(i+1)-1			% last but 1 pt on slice
				x1=sum(npt_tr(1:i))+j+2 ;
				x12=x1-npt_tr(i) ;
				x13=x12+1 ;
				x22=x13 ;
				x23=sum(npt_tr(1:i))+2+npt_qr ;
				ind_qr=[ind_qr [x1 x13 x12]' [x1 x23 x22]'] ;
			elseif j==npt_tr(i+1)		% last pt on slice
				x1=sum(npt_tr(1:i))+2+npt_qr ;
				x2=sum(npt_tr(1:(i-1)))+j+2 ;
				x3=x1-npt_tr(i) ;
				ind_qr=[ind_qr [x1 x3 x2]'] ;
			else            			% other pts on slice
				x1=sum(npt_tr(1:i))+j+2 ;
				x12=x1-npt_tr(i) ;
				x13=x12+1 ;
				x22=x13 ;
				x23=x1+1 ;
				ind_qr=[ind_qr [x1 x13 x12]' [x1 x23 x22]'] ;
			end
		end
	end
end


ntr_qr   =size(ind_qr,2) ; % nbr of triangle per 5th of sphere
[S_i,S_j]=find(ind_qr==1) ;
[B_i,B_j]=find(ind_qr==npt_sph) ;
[qs_i,qs_j]=find((ind_qr>(npt_qr+1))&(ind_qr<(3*npt_qr))) ;
ind_tr=[] ;

% shift all indices to cover the sphere
for i=0:4
	ind_tr = [ind_tr (ind_qr+i*npt_qr)] ;
	ind_tr(S_i,S_j+i*ntr_qr) = 1 ;
	ind_tr(B_i,B_j+i*ntr_qr) = npt_sph ;
end
for i=1:size(qs_i,1)
	ind_tr(qs_i(i),qs_j(i)+4*ntr_qr) = ...
		ind_tr(qs_i(i),(qs_j(i)+4*ntr_qr))-5*npt_qr ;
end
ntr_sph = size(ind_tr,2) ;% nbr of triangles on the sphere
	% or = 5/2*n^2
return
