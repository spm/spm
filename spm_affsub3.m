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
	free   = [ones(1,18) ones(1, 2)];

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
	free   = [ones(1,6) zeros(1,6) zeros(1, np)];

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
	free   = [ones(1,6) zeros(1,6) ones(1, np)];

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
	free   = [ones(1,6) zeros(1,6) ones(1, np)];

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
	free   = [ones(1,12) zeros(1, np)];

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
	free   = [ones(1,12) ones(1, np)];

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
	free   = [ones(1,12) ones(1, np)];
else
	error('I dont understand');
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
params = spm_affsub2(VG,VF, MG,MF, Hold,samp,params,free,pdesc,gorder);


% Tidy up
%-----------------------------------------------------------------------
for V = [VF VG]
	spm_unmap(V);
end
clear VF VG

