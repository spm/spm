function params = spm_affsub3(mode, PG, PF, Hold, samp, params,PW,MW)
% Highest level subroutine involved in affine transformations.
% FORMAT params = spm_affsub3(mode, PG, PF, Hold, samp, params)
%
% mode      - Mode of action.
% PG        - Matrix of template image name.
% PF        - Matrix of object image names.
% Hold      - Interpolation method.
% samp      - Frequency (in mm) of sampling.
% params    - Parameter estimates.
%
% optional:
% PW        - Name of weight image.
% MW        - Affine matrix mapping weight image to space of PG images.
%__________________________________________________________________________
%
% Currently mode must be one of the following:
%	'register1'
%		This is for use in Multimodal coregistration.
%		Each F is mapped to one G (with scaling), but the
%		rigid body components differ between the two sets
%		of registrations.
%	'rigid1'
%		Rigid body registration.
%		Each F is mapped to one G, without scaling.
%	'rigid2'
%		Rigid body registration.
%		Each F is mapped to one G, with scaling.
%	'rigid3'
%		Rigid body registration.
%		Each F is mapped to a linear combination of Gs.
%	'affine1'
%		Affine normalisation.
%		Each F is mapped to one G, without scaling.
%	'affine2'
%		Affine normalisation.
%		Each F is mapped to one G, with scaling.
%	'affine3'
%		Affine normalisation.
%		Each F is mapped to a linear combination of Gs.
%
%	'2d1'
%		For 2d rigid-body registration.
%__________________________________________________________________________
% %W% John Ashburner FIL %E%

if nargin<5 | nargin>8
	error('Incorrect usage.');
end

global sptl_Ornt

% Covariance matrix of affine normalization parameters. Mostly assumed
% to be diagonal, but covariances between the three zooms is also
% accounted for.
% Data determined from 51 normal brains.
%-----------------------------------------------------------------------
c11 = eye(3)*10000;                % Allow stdev of 100mm in all directions.
c22 = eye(3)*0.046;                % 7 degrees standard deviations.
c33 = [0.0021 0.0009 0.0013        % Covariances of zooms (from 51 brains).
       0.0009 0.0031 0.0014
       0.0013 0.0014 0.0024];
c44 = diag([0.18 0.11 1.79]*1e-3); % Variance of shears (from 51 brains).
pad = zeros(3);
covar = [c11 pad pad pad
         pad c22 pad pad
         pad pad c33 pad
         pad pad pad c44];
icovar = inv(covar);

ornt = sptl_Ornt;
ornt(7:9) = ornt(7:9).*[1.1 1.05 1.17];

nobayes = 0;

if strcmp(mode,'register1')
	% This is for use in Multimodal coregistration.
	% Each F is mapped to one G (with scaling), but the
	% rigid body components differ between the two sets
	% of registrations.
	%-----------------------------------------------------------------------
	pdesc   = [1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 0
		   0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 1]';
	gorder  = [1 2];
	free    = [ones(1,18) ones(1, 2)]';
	mean0   = [ornt([1:6 1:6 7:12]) 1 1]';
	icovar0 = zeros(length(mean0));
	icovar0( 1:6 , 1:6 ) = icovar( 1:6 , 1:6 );
	icovar0( 7:12, 7:12) = icovar( 1:6 , 1:6 );
	icovar0(13:18,13:18) = icovar( 7:12, 7:12);
	icovar0( 1:6 ,13:18) = icovar( 1:6 , 7:12);
	icovar0( 7:12,13:18) = icovar( 1:6 , 7:12);
	icovar0(13:18, 1:6 ) = icovar( 7:12, 1:6 );
	icovar0(13:18, 7:12) = icovar( 7:12, 1:6 );

elseif strcmp(mode,'rigid1')
	% Rigid body registration.
	% Each F is mapped to one G, without scaling.
	%-----------------------------------------------------------------------
	np     = size(PG,1);
	if np ~= size(PF,1)
		error('There should be the same number of object and template images');
	end

	pdesc   = [ones(12,np); eye(np)];
	gorder  = 1:np;
	free    = [ones(1,6) zeros(1,6) zeros(1, np)]';
	mean0   = [[0 0 0 0 0 0 1 1 1 0 0 0] ones(1,np)]';
	icovar0 = zeros(length(mean0));
	%icovar0(1:6,1:6)=icovar(1:6,1:6);
	nobayes = 1;

elseif strcmp(mode,'rigid2')
	% Rigid body registration.
	% Each F is mapped to one G, with scaling.
	%-----------------------------------------------------------------------
	np     = size(PG,1);
	if np ~= size(PF,1)
		error('There should be the same number of object and template images');
	end

	pdesc   = [ones(12,np); eye(np)];
	gorder  = 1:np;
	free    = [ones(1,6) zeros(1,6) ones(1, np)]';
	mean0   = [[0 0 0 0 0 0 1 1 1 0 0 0] ones(1,np)]';
	icovar0 = zeros(length(mean0));
	%icovar0(1:6,1:6)=icovar(1:6,1:6);

elseif strcmp(mode,'rigid3')
	% Rigid body registration.
	% Each F is mapped to a linear combination of Gs.
	%-----------------------------------------------------------------------
	np    = size(PG,1);
	if size(PF,1) ~= 1
		error('There should be one object image');
	end

	pdesc   = ones(12+np,1);
	gorder  = ones(1,np);
	free    = [ones(1,6) zeros(1,6) ones(1, np)]';
	mean0   = [[0 0 0 0 0 0 1 1 1 0 0 0] ones(1,np)]';
	icovar0 = zeros(length(mean0));
	%icovar0(1:6,1:6)=icovar(1:6,1:6);
	nobayes = 1;

elseif strcmp(mode,'2d1')
	% Rigid body registration.
	% Each F is mapped to a linear combination of Gs.
	%-----------------------------------------------------------------------
	np    = size(PG,1);
	if size(PF,1) ~= 1
		error('There should be one object image');
	end
	pdesc   = ones(12+np,1);
	gorder  = ones(1,np);
	free    = [[1 1 0 0 0 1] zeros(1,6) ones(1, np)]';
	mean0   = [[0 0 0 0 0 0 1 1 1 0 0 0] ones(1,np)]';
	icovar0 = zeros(length(mean0));
	%icovar0(1:6,1:6)=icovar(1:6,1:6);
	nobayes = 1;

elseif strcmp(mode,'affine1')
	% Affine normalisation.
	% Each F is mapped to one G, without scaling.
	%-----------------------------------------------------------------------
	np     = size(PG,1);
	if np ~= size(PF,1)
		error('There should be the same number of object and template images');
	end

	pdesc   = [ones(12,np); eye(np)];
	gorder  = 1:np;
	free    = [ones(1,12) zeros(1, np)]';
	mean0   = [ornt ones(1,np)]';
	icovar0 = zeros(length(mean0));
	icovar0(1:12,1:12) = icovar(1:12,1:12);

elseif strcmp(mode,'affine2')
	% Affine normalisation.
	% Each F is mapped to one G, with scaling.
	%-----------------------------------------------------------------------
	np     = size(PG,1);
	if np ~= size(PF,1)
		error('There should be the same number of object and template images');
	end

	pdesc   = [ones(12,np); eye(np)];
	gorder  = 1:np;
	free    = [ones(1,12) ones(1, np)]';
	mean0   = [ornt ones(1,np)]';
	icovar0 = zeros(length(mean0));
	icovar0(1:12,1:12) = icovar(1:12,1:12);

elseif strcmp(mode,'affine3')
	% Affine normalisation.
	% Each F is mapped to a linear combination of Gs.
	%-----------------------------------------------------------------------
	np    = size(PG,1);
	if size(PF,1) ~= 1
		error('There should be one object image');
	end

	pdesc   = ones(12+np,1);
	gorder  = ones(1,np);
	free    = [ones(1,12) ones(1, np)]';
	mean0   = [ornt ones(1,np)]';
	icovar0 = zeros(length(mean0));
	icovar0(1:12,1:12) = icovar(1:12,1:12);

else
	error('I dont understand');
end


% Map the images & get their positions in space.
%-----------------------------------------------------------------------
VG = spm_vol(PG);
VF = spm_vol(PF);
if nargin<8,
	VW = [];
else,
	if ~isempty(PW),
		VW = spm_vol(PW);
		for i=1:length(VW)
			VW(i).mat = MW*VW(i).mat;
		end;
	else,
		VW = [];
	end;
end

% Do the optimisation
%-----------------------------------------------------------------------
if nargin == 5
	params = mean0;
end

if nobayes == 1
	[params] = spm_affsub2(VG,VF,VW, Hold,samp,params,free,pdesc,gorder);
else
	[params] = spm_affsub2(VG,VF,VW, Hold,samp,params,free,pdesc,gorder,mean0,icovar0);
end

