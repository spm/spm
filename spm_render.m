function spm_render(XYZ,t,V)
% Render blobs on surface of a 'standard' brain
% FORMAT spm_render
%_______________________________________________________________________
% 
% FORMAT spm_render(XYZ,t,V)
% XYZ - the x, y & z coordinated of the transformed t values.
% t   - the transformed t values
% V   - a vector containing the voxel sizes in positions 5-7.
%_______________________________________________________________________
% 
% spm_render prompts for details of a SPM{Z} that is then displayed
% superimposed on the surface of a standard brain.
% The blobs which are displayed are the integral of all transformed t
% values which are less than 3cm deep.


% %W% John Ashburner FIL %E%


if (nargin==0)
	% User interface
	%_______________________________________________________________________
	spm_figure('Clear','Interactive');
	tmp  = spm_get(1,'.mat','select SPMt.mat for analysis','SPMt');
	d    = strrep(tmp,'/SPMt.mat','');
	K    = [];

	set(spm_figure('FindWin','Interactive'),'Name','loading','Pointer','watch');
	load([d,'/SPM' ])
	load([d,'/XYZ' ])
	load([d,'/SPMt'])
	spm_figure('Clear','Interactive');

	con  = 0;
	while any(con < 1 | con > size(CONTRAST,1))
		con = spm_input(sprintf('contrast[s] ? 1 - %i',size(CONTRAST,1)),1);
	end

	% get height threshold [default = 3.2]
	%---------------------------------------------------------------------------
	U    = spm_input('height threshold {Z value}',2,'e',3.2);

	% get extent threshold [default = E{n} - expected voxels per cluster]
	% Omit spatial extent threshold for multiple contrasts.
	%---------------------------------------------------------------------------
	if length(con) == 1
		[P,EN,Em,En,Pk] = spm_P(1,W,U,0,S);
		k    = spm_input('extent threshold {voxels}',3,'e',round(En));
	else
		k    = 0;
	end


	% accommodate masking if required
	%---------------------------------------------------------------------------
	if length(con) > 1
		Q = all(SPMt(con,:) > U);

		c = CONTRAST(con,:);
		g = [K H C B G];
		g = c*pinv(g'*g)*c';
		r = inv(diag(sqrt(diag(g))))'*g*inv(diag(sqrt(diag(g))));
		t = sum(SPMt(con,Q))/sqrt(sum(r(:)));
	else
		Q = SPMt(con,:) > U;
		t = SPMt(con,Q);
	end

	% return if there are no voxels
	%---------------------------------------------------------------------------
	if sum(Q) == 0
		figure(spm_figure('FindWin','Graphics'));
		spm_figure('Clear','Graphics');
		axis off
		text(0,0.8,d);
		text(0,0.7,'No voxels above this threshold {u}','FontSize',16);
		return
	end

	XYZ  = XYZ(:,Q);


	% apply threshold {k}
	%---------------------------------------------------------------------------
	A         = spm_clusters(XYZ,V([4 5 6]));
	Q         = [];
	for i     = 1:max(A)
		j = find(A == i);
		if length(j) >= k
			Q = [Q j];
		end
	end

	% return if there are no voxels
	%---------------------------------------------------------------------------
	if sum(Q) == 0
		figure(spm_figure('FindWin','Graphics'));
		spm_figure('Clear','Graphics');
		axis off
		text(0,0.8,d);
		text(0,0.7,'No clusters above this threshold {k}','FontSize',16);
		return
	end

	t    = t(Q);
	XYZ  = XYZ(:,Q);

	spm_render(XYZ,t,V);
	return;
end



% Perform the rendering
%_______________________________________________________________________

set(spm_figure('FindWin','Interactive'),'Name','executing','Pointer','watch');


if (nargin<3)
	V = [0 0 0 2 2 2];
end

load render
load Split
colormap(split);

mx = 0;

spm_progress_bar('Init', 8, 'Making pictures', 'Number completed');

for i=1:size(Matrixes,1)
	eval(['MM = ' Matrixes(i,:) ...
		'; ren = ' Rens(i,:) '; dep = ' Depths(i,:) ';']);

	% transform from Taliarach space to space of the rendered image
	%---------------------------------------------------------------------------
	xyz      = (MM(1:3,1:3)*XYZ);
	xyz(1,:) = xyz(1,:) + MM(1,4);
	xyz(2,:) = xyz(2,:) + MM(2,4);
	xyz(3,:) = xyz(3,:) + MM(3,4);
	xyz      = round(xyz);

	% field should be sampled every mm.
	%---------------------------------------------------------------------------
	t1       = [];
	xyz1     = [];
	xx       = xyz(1,:);
	for x1=(-V(4)/2):(V(4)/2-eps*20)
		xyz(1,:) = xx+x1;
		xyz1 = [xyz1 xyz];
		t1   = [t1 t];
	end
	t0        = [];
	yy        = xyz1(2,:);
	xyz       = [];
	for y1=(-V(5)/2):(V(5)/2-eps*20)
		xyz1(2,:) = yy+y1;
		xyz = [xyz xyz1];
		t0  = [t0 t1];
	end


	% only use values which will fit on the image
	%---------------------------------------------------------------------------
	msk = find((xyz(1,:) >= 1) & (xyz(1,:) <= size(dep,1)) ...
		&  (xyz(2,:) >= 1) & (xyz(2,:) <= size(dep,2)));
	xyz      = xyz(:,msk);
	t0       = t0(msk);

	X = zeros(size(dep));

	if ~isempty(xyz)

		% calculate 'depth' of values, and ignore any which are more than 30 mm
		% behind the surface.
		%---------------------------------------------------------------------------
		z1  = dep(xyz(1,:)+(xyz(2,:)-1)*size(dep,1));
		msk = find(xyz(3,:) < (z1+30) & xyz(3,:) > (z1-5));
		xyz = xyz(:,msk);
		t0  = t0(msk);

		% generate an image of the integral of the blob values.
		%---------------------------------------------------------------------------
		if ~isempty(xyz)
			X = full(sparse(xyz(1,:), xyz(2,:), t0, size(dep,1), size(dep,2)));
		end
	end

	mx = max([mx max(max(X))]);
	eval(['X' num2str(i) ' = X;']);
	spm_progress_bar('Set', i);
end

spm_progress_bar('Clear');
figure(spm_figure('FindWin','Graphics'));
spm_figure('Clear','Graphics');

% Combine the brain surface renderings with the blobs, and display.
%---------------------------------------------------------------------------
for i=1:size(Matrixes,1)
	eval(['ren = ' Rens(i,:) '; X = X' num2str(i) ';']);
	msk = find(X ~= 0);
	ren(msk) = X(msk)*(64/(mx))+64;
	subplot(4,2,i);
	image(ren);
	axis image off xy
end

spm_figure('Clear','Interactive');

