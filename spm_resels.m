function [R] = spm_resels(FWHM,L,SPACE)
% Returns the RESEL counts of a search volume
% FORMAT [R] = spm_resels(FWHM,L,SPACE)
% FWHM       - smoothness of the component fields {FWHM - voxels}
% L          - space definition            {in voxels}
%                L = radius                {Sphere}
%                L = [height width length] {Box}
%                L = XYZ pointlist         {Discrete voxels}
%                L = Mapped image volume   {Image}
% SPACE      - Search space
%               'S' - Sphere
%               'B' - Box
%               'V' - Discrete voxels
%               'I' - Image VOI
%
% R          - RESEL counts {adimensional}
%
%___________________________________________________________________________
% For one or two dimensional spaces the appropriate manifold is
% used (e.g. sphere -> disc -> line).  
%
% Reference : Worsley KJ et al 1996, Hum Brain Mapp. 4:58-73
%
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston & Matthew Brett
% $Id: spm_resels.m 112 2005-05-04 18:20:52Z john $


% Dimensionality
%---------------------------------------------------------------------------
switch SPACE

case 'S'                                                            % Sphere
	%-------------------------------------------------------------------
	s     = L(:)./FWHM(:);
        s     = s(s > 0);
        if length(s) == 2,  SPACE = 'D';  end
        if length(s) == 1,  SPACE = 'L';  end

case 'B'                                                               % Box
	%-------------------------------------------------------------------
	s     = L(:)./FWHM(:);
        s     = s(s > 0);
        if length(s) == 2,  SPACE = 'R';  end
        if length(s) == 1,  SPACE = 'L';  end
end


% Default {sphere - assuming L = volume i.e. number of voxels)
%---------------------------------------------------------------------------
if nargin < 3
	SPACE = 'S';
	L     = (L*(3/4)/pi)^(1/3);
end


% RESEL Counts (R)
%===========================================================================

switch SPACE

case 'S'                                                            % Sphere
	%-------------------------------------------------------------------
	s     = prod(s).^(1/3);
	R     = [1 4*s 2*pi*s^2 (4/3)*pi*s^3];

case 'D'                                                              % Disc
	%-------------------------------------------------------------------
	s     = prod(s).^(1/2);
	R     = [1 pi*s pi*s^2 0];

case 'B'                                                               % Box
	%-------------------------------------------------------------------
	R     = [1 sum(s) (s(1)*s(2) + s(2)*s(3) + s(1)*s(3)) prod(s)];

case 'R'                                                         % Rectangle
	%-------------------------------------------------------------------
	R     = [1 sum(s) prod(s) 0];

case 'L'                                                              % Line
	%-------------------------------------------------------------------
	R     = [1 s 0 0];

case 'V'                                                            % Voxels
	%-------------------------------------------------------------------
	R     = spm_Pec_resels(L,FWHM);

case 'I'                                                             % Image
	%-------------------------------------------------------------------
	R     = spm_resels_vol(L,FWHM);
	R     = R' .* [1 2/3 2/3 1];	       %-KJW "knobliness" correction
end

