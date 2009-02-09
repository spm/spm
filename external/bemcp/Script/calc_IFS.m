% function calc_IFS

% function calc_IFS
%
% Calculates the IFS (Intermediate Forward Solution) matrices
% using the results of prep_C function.
% 
% Loads Cp1, Cp2 and param
% Saves IFS

TICTOC 		= clock ;
TICTOCcpu	= cputime;
save meas_time

normal = 1 ;

disp('flag 3.1.1')
load C11_21
if normal
	tmp1 = C21/C11st ;
else
	n = size(C11st,1)
	n2 = floor(n/2)
	CAi = C11st((n2+1):end,1:n2)/C11st(1:n2,1:n2) ;
	F = C11st((n2+1):end,(n2+1):end) - CAi*C11st(1:n2,(n2+1):end) ;
	UAi = C21(:,1:n2)/C11st(1:n2,1:n2) ;
	VFi = C21(:,(n2+1):end)/F ;
	clear C21
	UAiBFi = UAi*(C11st(1:n2,(n2+1):end)/F) ;
	clear C11st
	tmp1 = [((UAiBFi*CAi)+UAi)-(VFi*CAi) -UAiBFi+VFi] ;
end

save tmp1 tmp1
clear all


disp('flag 3.1.2')
load C12_22
tmp2 = C12/C22st ;
save tmp2 tmp2
clear tmp2

disp('flag 3.1.3')
load C32
load tmp1
tmp3 = C32/((-tmp1*C12) + C22st ) ;
save tmp3
clear all

disp('flag 3.1.4')
load C11_21
load tmp2
if normal
	load C31
	tmp4 = C31/((-tmp2*C21) + C11st ) ;
else
	C11tmp = -tmp2*C21 ;
	save C11tmp C11tmp
	save C11_21 C11st
	clear all
	load C11tmp
	load C11_21
	C11bis = C11tmp + C11st ;
	save C11_21 C11bis
	clear all
	delete C11tmp.mat
	load C11_21
	load C31
	n = size(C11bis,1)
	n2 = floor(n/2)
	CAi = C11bis((n2+1):end,1:n2)/C11bis(1:n2,1:n2) ;
	F = C11bis((n2+1):end,(n2+1):end) - CAi*C11bis(1:n2,(n2+1):end) ;
	UAi = C31(:,1:n2)/C11bis(1:n2,1:n2) ;
	VFi = C31(:,(n2+1):end)/F ;
	clear C31
	UAiBFi = UAi*(C11bis(1:n2,(n2+1):end)/F) ;
	clear C11bis
	tmp4 = [((UAiBFi*CAi)+UAi)-(VFi*CAi) -UAiBFi+VFi] ;
end
save tmp4 tmp4
clear all

delete C11_21.mat C12_22.mat C32.mat C31.mat

disp('flag 3.1.5')
load tmp1
load tmp3
load tmp4
tmp5 = (tmp3*tmp1)-tmp4 ;
save tmp5 tmp5
clear all

delete tmp1.mat

disp('flag 3.1.6')
load tmp2
load tmp3
load tmp4
tmp6 = (tmp4*tmp2)-tmp3 ;
save tmp6 tmp6
clear all

delete tmp2.mat tmp3.mat tmp4.mat

disp('flag 5.1')
load C13_23_33
load tmp5
load tmp6
gama1 = - C33st - (tmp5*C13) - (tmp6*C23) ;
save gama1 gama1
clear C33st C23 C13
delete C13_23_33.mat 

load Param

A = gama1(1:nRem_vert,1:nRem_vert) ;
B = gama1(1:nRem_vert,nRem_vert+1:tsph3.nr(1)) ;
C = gama1(nRem_vert+1:tsph3.nr(1),1:nRem_vert) ;
D = gama1(nRem_vert+1:tsph3.nr(1),nRem_vert+1:tsph3.nr(1)) ;
clear gama1
delete gama1.mat

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
save IFS IFS1 IFS2 IFS3
clear D B tmp5 tmp6 gama_el
delete tmp5.mat tmp6.mat

load meas_time
elapsed_time = etime(clock, TICTOC)
elapsed_cpu_time = cputime-TICTOCcpu
delete meas_time.mat


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load C11_21
%tmp1 = C21/C11st ;
n = size(C11st,1)
n2 = floor(n/2)

CAi = C11st((n2+1):end,1:n2)/C11st(1:n2,1:n2) ;
F = C11st((n2+1):end,(n2+1):end) - CAi*C11st(1:n2,(n2+1):end) ;
UAi = C21(:,1:n2)/C11st(1:n2,1:n2) ;
VFi = C21(:,(n2+1):end)/F ;
clear C21
UAiBFi = UAi*(C11st(1:n2,(n2+1):end)/F) ;
clear C11st
tmp1 = [((UAiBFi*CAi)+UAi)-(VFi*CAi) -UAiBFi+VFi] ;


