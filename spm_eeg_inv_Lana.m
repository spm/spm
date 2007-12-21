function [Lan,nit]=spm_eeg_inv_Lana(XYZmm,SseXYZ,Rsc,Rsk,Rbr,sigma)

% FUNCTION Lan,nit]=spm_eeg_inv_Lana(XYZmm,SseXYZ,Rsc,Rsk,Rbr)
% Calculates the leadfield in a 3-shell sphere head model for a set of
% distributed dipoles. As this is a spherical model, the solution is analytical
% but it relies on a truncated infinite series.
% IN :
%   - XYZmm         : dipoles location, moved to fit the spherical model.
%   - SseXYZ        : coordinates of the electrodes on the "scalp sphere".
%   - Rsc,Rsk,Rbr   : radii of the scalp, skull, brain spheres.
%   - sigma         : conductivities for the scalp, skull and brain volumes.
% OUT :
%   - Lan : (orientation free) leadfield obtained with the analytical formula
%   - nit : number of terms used in the infinite series
%
% Based on formulas found in :
% James P. Ary, Stanley A. Klein, Derek H. Fender
% Location of fources of evoked scalp potentials: Corrections for skull
% and scalp thicknesses
% IEEE TBME, 28-6, 447452, 1981
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips,
% $Id: spm_eeg_inv_Lana.m 1039 2007-12-21 20:20:38Z karl $

if nargin<6
    sigma = [.33 .004 .33];
end
e=sigma(2)/sigma(1) ;
f1=Rbr/Rsc ;
f2=Rsk/Rsc ;

Ndip = size(XYZmm,2) ;
Lan=zeros(size(SseXYZ,2),3*Ndip) ;
nit =zeros(1,Ndip) ;

% aff=1000 ;
% fprintf('\n')
for ii=1:Ndip
	%	disp(['ii = ',num2str(ii)])
% 	if rem(ii,aff)==0
%         fprintf(' %3.1f%% calculated. \n',ii/Ndip*100)
%     end
	xyz_d = XYZmm(:,ii) ;
    % This is not the same spherical coordinates as produced by cart2sph !!!
    % azimuth is expressed as:
    % - in here: phi = angle betwee e_x and proj of point on 'oxy' plane
    % - in cart2sph: *th* = idem
    % elevation is expressed as:
    % - in here: thet = angle between e_z and point 
    % - in cartsph : *phi* = angle between point and 'oxy' plane
    %                i.e. *phi* = pi/2 - thet
    
	no = norm(xyz_d) ;
	b=no/Rsc ; % eccentricity
	if no==0
		thet=0 ;
	else
		thet=acos(xyz_d(3)/no) ;
	end
	phi = atan2(xyz_d(2),xyz_d(1)) ;

	% Rotation matrix in order to bring the point (i.e. dipole) on e_z
    rotyz = [cos(phi)*cos(thet) sin(phi)*cos(thet) -sin(thet) ;...
		     -sin(phi)          cos(phi)            0 ;...
		     cos(phi)*sin(thet) sin(phi)*sin(thet) cos(thet)] ;

    % Apply the same rotation on electrode locations
	Pt_m = rotyz * SseXYZ  ;
    % and get their spherical coordinates.
	alpha_p = acos(Pt_m(3,:)/Rsc)' ;
	beta_p = atan2(Pt_m(2,:),Pt_m(1,:))' ;

	mx = [ cos(phi)*cos(thet) -sin(phi) sin(thet)*cos(phi) ] ;
	my = [ sin(phi)*cos(thet) cos(phi) sin(thet)*sin(phi)] ;
	mz = [ -sin(thet) 0 cos(thet) ] ;
    
    % All data are prepared proceed to main calculation for all electrodes at once.
	[VR,VT,nit(ii)] = V1dip(Rsc,f1,f2,e,sigma(1),b,alpha_p) ;
    y  = cos(beta_p) ; yy = sin(beta_p) ;
    
    % the source orientation (after rotation) m_xyz is combined here
    % with the potential V for tangent (VT) and radial (VR) dipole
    % and electrodes location beta_p).
	Lan(:,3*ii-2) = mx(1)*y.*VT + mx(2)*yy.*VT + mx(3)*VR ;
	Lan(:,3*ii-1) = my(1)*y.*VT + my(2)*yy.*VT + my(3)*VR ;
	Lan(:,3*ii) = mz(1)*y.*VT + mz(3)*VR ;
end
% fprintf('\n')
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SUBFUNCTION
%
function [VR,VT,n] = V1dip(R,f1,f2,e,sig1,b,alpha_p)

x = cos(alpha_p) ; xx = sin(alpha_p) ; 
lim=1e-12 ; nlim=200 ;

% Formula (2a) for d_n
d1 = e*e + 1 + 2.5*e + (1-2*e*e+e)*(f1^3-f2^3) - (1-e)^2*(f1/f2)^3 ;     % n=1
d2 = 2*e*e + 2 + e*13/3 + (2-3*e*e+e)*(f1^5-f2^5) - 2*(1-e)^2*(f1/f2)^5 ;% n=2

% This correspnonds to the terms using 'n', 'b' & 'e' in (2)
% (2n+1)/n * e*(2n+1)^2 * b^(n-1) /(d_n*(n+1)) = (2n+1)^3 / ( n*d_n*(n+1) ) * e/b^(n-1)
w1 = 27/2*e/d1 ;    % n=1
w2 = 125/6*b*e/d2 ; % n=2

P0_1 = x ; P1_1 = xx ;
P0_2 = .5*(3*x.*x-1) ; P1_2 = 3*xx.*x ;

VtempR = w1*P0_1 + w2*2*P0_2 ;
VtempT = w1*P1_1 + w2*P1_2 ;

Pn_2=P0_1 ; Pn_1=P0_2 ;
P1n_2=P1_1 ; P1n_1=P1_2 ;
VtempR_1=VtempR ;
VtempT_1=VtempT ;

n=3 ; diR=10 ; diT=10 ;

while ((max(abs(diR))>lim) & (max(abs(diT))>lim) & (n<nlim))
    % Update Legendre polynome.
	Pn  = (2*n-1)/n*x.*Pn_1-(n-1)/n*Pn_2 ;
	P1n = -(2*n-1)/(1-n)*x.*P1n_1+n/(1-n)*P1n_2 ;
    % Equation (2a) as in paper
	dn =  ( (n+1)*e + n )*( (n*e)/(n+1) + 1 ) ...
		+ (1-e)*( (n+1)*e + n ) * ( f1^(2*n+1) - f2^(2*n+1) )...
		- n*(1-e)^2*(f1/f2)^(2*n+1) ;
	wn = b^(n-1) * e * (2*n+1)^3 / dn / (n+1) / n ;

	VtempR = VtempR_1 + wn*n*Pn ;
	VtempT = VtempT_1 + wn*P1n ;

	diR=(VtempR-VtempR_1) ; % ./ VtempR_1 ;
	diT=(VtempT-VtempT_1) ; % ./ VtempT_1 ;

	n=n+1 ;
	Pn_2 = Pn_1 ; P1n_2 = P1n_1 ;
	Pn_1 = Pn   ; P1n_1 = P1n ;
	VtempR_1 = VtempR ;
	VtempT_1 = VtempT ;
end

VR = VtempR / (4*pi*sig1*R^2) ;
VT = VtempT / (4*pi*sig1*R^2) ;

return