function params = spm_affsub3(mode, PG, PF, Hold, samp, params)
% Highest level subroutine involved in affine transformations.
% FORMAT params = spm_affsub3(mode, PG, PF, Hold, samp, params)
%
% mode      - Mode of action.
% PG        - Matrix of template image name.
% PF        - Matrix of object image names.
% Hold      - Interpolation method.
% samp      - Frequency (in mm) of sampling.
% params    - Parameter estimates.
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
%		Each F is mapped to one G, without (much) scaling.
%	'rigid2'
%		Rigid body registration.
%		Each F is mapped to one G, with scaling.
%	'rigid3'
%		Rigid body registration.
%		Each F is mapped to a linear combination of Gs.
%	'affine1'
%		Affine normalisation.
%		Each F is mapped to one G, without (much) scaling.
%	'affine2'
%		Affine normalisation.
%		Each F is mapped to one G, with scaling.
%	'affine3'
%		Affine normalisation.
%		Each F is mapped to a linear combination of Gs.
%
%__________________________________________________________________________
% %W% John Ashburner FIL %E%

if nargin<5 | nargin>6
	error('Incorrect usage.');
end

if strcmp(mode,'register1')
	% This is for use in Multimodal coregistration.
	% Each F is mapped to one G (with scaling), but the
	% rigid body components differ between the two sets
	% of registrations.
	%-----------------------------------------------------------------------
	pdesc  = [1 1 1 1 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 0
		  0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 0 1]';
	gorder = [1 2];
	covar = diag([...
		1000 1000 1000 30*pi/180 10*pi/180 20*pi/180 ...
		1000 1000 1000 30*pi/180 10*pi/180 20*pi/180 ...
		0.2 0.2 0.2 0.02 0.02 0.02 ...
		1e12 1e12].^2);
	mean0  = [0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 1 1]';

elseif strcmp(mode,'rigid1')
	% Rigid body registration.
	% Each F is mapped to one G, without (much) scaling.
	%-----------------------------------------------------------------------
	np     = size(PG,1);
	if np ~= size(PF,1)
		error('There should be the same number of object and template images');
	end

	pdesc  = [ones(12,np); eye(np)];
	gorder = 1:np;
	covar  = diag([[1000 1000 1000 30*pi/180 10*pi/180 20*pi/180 ...
		eps eps eps eps eps eps] ones(1, np)*0.05].^2);
	mean0  = [[0 0 0 0 0 0 1 1 1 0 0 0] ones(1,np)]';

elseif strcmp(mode,'rigid2')
	% Rigid body registration.
	% Each F is mapped to one G, with scaling.
	%-----------------------------------------------------------------------
	np     = size(PG,1);
	if np ~= size(PF,1)
		error('There should be the same number of object and template images');
	end

	pdesc  = [ones(12,np); eye(np)];
	gorder = 1:np;
	covar  = diag([[1000 1000 1000 30*pi/180 10*pi/180 20*pi/180 ...
		eps eps eps eps eps eps] ones(1, np)*1e12].^2);
	mean0  = [[0 0 0 0 0 0 1 1 1 0 0 0] ones(1,np)]';

elseif strcmp(mode,'rigid3')
	% Rigid body registration.
	% Each F is mapped to a linear combination of Gs.
	%-----------------------------------------------------------------------
	np    = size(PG,1);
	if size(PF,1) ~= 1
		error('There should be one object image');
	end

	pdesc  = ones(12+np,1);
	gorder = ones(1,np);
	covar  = diag([[1000 1000 1000 30*pi/180 10*pi/180 20*pi/180 ...
		eps eps eps eps eps eps] ones(1, np)*1e12].^2);
	mean0  = [[0 0 0 0 0 0 1 1 1 0 0 0] ones(1,np)]';

elseif strcmp(mode,'affine1')
	% Affine normalisation.
	% Each F is mapped to one G, without (much) scaling.
	%-----------------------------------------------------------------------
	np     = size(PG,1);
	if np ~= size(PF,1)
		error('There should be the same number of object and template images');
	end

	pdesc  = [ones(12,np); eye(np)];
	gorder = 1:np;
	covar  = diag([[1000 1000 1000 30*pi/180 10*pi/180 20*pi/180 ...
		0.2 0.2 0.2 0.02 0.02 0.02] ones(1, np)*0.05].^2);
	mean0  = [[0 0 0 0 0 0 1 1 1 0 0 0] ones(1,np)]';

elseif strcmp(mode,'affine2')
	% Affine normalisation.
	% Each F is mapped to one G, with scaling.
	%-----------------------------------------------------------------------
	np     = size(PG,1);
	if np ~= size(PF,1)
		error('There should be the same number of object and template images');
	end

	pdesc  = [ones(12,np); eye(np)];
	gorder = 1:np;
	covar  = diag([[1000 1000 1000 30*pi/180 30*pi/180 30*pi/180 ...
		0.2 0.2 0.2 0.02 0.02 0.02] ones(1, np)*1e12].^2);
	mean0  = [[0 0 0 0 0 0 1 1 1 0 0 0] ones(1,np)]';

elseif strcmp(mode,'affine3')
	% Affine normalisation.
	% Each F is mapped to a linear combination of Gs.
	%-----------------------------------------------------------------------
	np    = size(PG,1);
	if size(PF,1) ~= 1
		error('There should be one object image');
	end

	pdesc  = ones(12+np,1);
	gorder = ones(1,np);
	covar  = diag([[1000 1000 1000 30*pi/180 10*pi/180 20*pi/180 ...
		0.2 0.2 0.2 0.02 0.02 0.02] ones(1, np)*1e12].^2);
	mean0  = [[0 0 0 0 0 0 1 1 1 0 0 0] ones(1,np)]';
else
	error('I dont understand');
end

if nargin==5
	params = mean0;
end

% Map the images & get their positions in space.
%-----------------------------------------------------------------------
for i=1:size(PG,1)
	filename = deblank(PG(i,:));
	VG(:,i) = spm_map(filename);
	MG(:,i) = reshape(spm_get_space(filename),16,1);
end
for i=1:size(PF,1)
	filename = deblank(PF(i,:));
	VF(:,i) = spm_map(filename);
	MF(:,i) = reshape(spm_get_space(filename),16,1);
end


% Do the optimisation
%-----------------------------------------------------------------------
params = spm_affsub2(VG,VF, MG,MF, Hold,samp,params,mean0,covar,pdesc,gorder);


% Tidy up
%-----------------------------------------------------------------------
for V = [VF VG]
	spm_unmap(V);
end
clear VF VG

