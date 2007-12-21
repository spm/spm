function varargout = spm_eeg_inv_vert2Rsph(flag,spheres,varargin)

% Transform coordinates (in mm) of any point into the coordinates 
% for a "Realistic sphere". (flag = 1)
%
% or does it the other way round, i.e. 
%
% Transforms the coordintates of a point in the "Realistic sphere" space
% back to the real world space (in mm). (flag = 2)
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips,
% $Id: spm_eeg_inv_vert2Rsph.m 1039 2007-12-21 20:20:38Z karl $

if flag==1
    vert = varargin{1};

	Nvert = size(vert,2);
	[th,phi,r_d] = cart2sph(vert(1,:)-spheres.centre(1), ...
                                vert(2,:)-spheres.centre(2), ...
                                vert(3,:)-spheres.centre(3));
	vTPR  = [th ; phi ; r_d];
	[XSv,YSv,ZSv] = sph2cart(vTPR(1,:), vTPR(2,:), 1);
	SvXYZ = [XSv;YSv;ZSv] ;

	% Interpolation of the scalp radius for each vert,
	% throug a k-nearest neighbour approach.
	Nneighb = 4;
	m = -2;
	% 'radius' of scalp in direction of each dipole
	r_ds = zeros(Nvert,1);
	for ii=1:Nvert
        [cos_angl,perm] = sort(abs(acos(SvXYZ(:,ii)'*spheres.Ss_XYZ)));
        sc_i  = perm(1:Nneighb);
        d_neighb = sqrt(sum((SvXYZ(:,ii)*ones(1,Nneighb)-spheres.Ss_XYZ(:,sc_i)).^2));
        r_ds(ii) = sum(d_neighb.^m .* spheres.r_s(sc_i))/sum(d_neighb.^m);
	end
	
	Svert = SvXYZ*diag(r_d'.*(spheres.R./r_ds));
    varargout{1} = Svert;
    varargout{2} = r_ds;
elseif flag==2
     Svert = varargin{1};
	Nvert = size(Svert,2);
	r_Svert = sqrt(sum(Svert.^2));
	sSvert = Svert*diag(1./r_Svert);

	% Interpolation of the scalp radius for each vert,
	% throug a k-nearest neighbour approach.
	Nneighb = 4;
	m = -2;
	% 'radius' of scalp in direction of each dipole
	r_ds = zeros(Nvert,1);
	for ii=1:Nvert
        [cos_angl,perm] = sort(abs(acos(sSvert(:,ii)'*spheres.Ss_XYZ)));
        sc_i  = perm(1:Nneighb);
        d_neighb = sqrt(sum((sSvert(:,ii)*ones(1,Nneighb)-spheres.Ss_XYZ(:,sc_i)).^2));
        r_ds(ii) = sum(d_neighb.^m .* spheres.r_s(sc_i))/sum(d_neighb.^m);
	end
	
	Rr_ds = r_Svert.*(r_ds'/spheres.R);
	vert = sSvert*diag(Rr_ds);
	vert = vert+spheres.centre*ones(1,Nvert);
    varargout{1} = vert;
    varargout{2} = r_ds;
else
    error('Wrong flag specified. Don''t know what to do...');
end

