function [R] = spm_resels(FWHM,L,SPACE)
% Returns the RESEL counts of a search volume
% FORMAT [R] = spm_resels(FWHM,L,SPACE)
% FWHM       - smoothness of the component fields {FWHM - voxels}
% L          - space definition            {in voxels}
%                L = radius                {Sphere}
%                L = [height width length] {Box}
%                L = XYZ pointlist         {Discrete voxels}
% SPACE      - Search space
%               'S' - Sphere
%               'B' - Box
%               'V' - Discrete voxels
%
% R          - RESEL counts {adimensional}
%
%___________________________________________________________________________
%
% Reference : Worsley KJ et al 1996, Hum Brain Mapp. 4:58-73
%
%___________________________________________________________________________
% %W% Karl Friston %E%

% Dimensionality
%---------------------------------------------------------------------------
D     = length(FWHM);

% Default {sphere - assuming L = volume)
%---------------------------------------------------------------------------
if nargin < 3
	SPACE = 'S';
	L     = (L*(3/4)/pi)^(1/3);
end


% RESEL Counts (R)
%===========================================================================
if      SPACE == 'S'

	% Sphere
	%-------------------------------------------------------------------
	s     = L(:)./FWHM(:);
	s     = prod(s).^(1/D);
	R     = [1 4*s 2*pi*s^2 (4/3)*pi*s^3];

elseif  SPACE == 'B'

	% Box
	%-------------------------------------------------------------------
	s     = L(:)./FWHM(:);
	R     = [1 sum(s) (s(1)*s(2) + s(2)*s(3) + s(1)*s(3)) prod(s)];

elseif  SPACE == 'V'

	% Voxels
	%-------------------------------------------------------------------
	R     = spm_Pec_resels(L,FWHM);

end

