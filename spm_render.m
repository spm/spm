function spm_render(dat)
% Render blobs on surface of a 'standard' brain
% FORMAT spm_render
%_______________________________________________________________________
% 
% FORMAT spm_render(dat)
% dat - a vertical cell array of length 1 to 3
%       - each element is a structure containing:
%         - XYZ - the x, y & z coordinated of the transformed t values.
%         - t   - the transformed t values
%         - v   - a vector containing the voxel sizes.
%
% Without arguments, spm_render acts as its own UI.
%_______________________________________________________________________
% 
% spm_render prompts for details of up to three SPM{Z}s or SPM{t}s that
% are then displayed superimposed on the surface of a standard brain.
% The first is shown in red, then green then blue.
% The blobs which are displayed are the integral of all transformed t
% values, exponentially decayed according to their depth. Voxels that
% are 10cm behind the surface have half the intensity of ones at the
% surface.
%_______________________________________________________________________
% %W% John Ashburner FIL %E%


if (nargin==0)

	spm_figure('Clear','Interactive');
	Finter = spm_figure('FindWin','Interactive');
	set(Finter,'Name','SPM render')
	global CWD

	% Which SPM
	%-----------------------------------------------------------------------
	SPMZ     = spm_input('which SPM',1,'b','SPM{Z}|SPM{F}',[1 0]);
	SPMF     = ~SPMZ;

	num = spm_input('Number of sets',2,'1 set|2 sets|3 sets',[1 2 3]);

	dat = cell(num,1);
	for i=1:num
		
		% Get thresholded data, thresholds and parameters
		%-----------------------------------------------------------------------
		if SPMZ
			[t,XYZ,QQ,U,k,s,w] = spm_projections_ui('Results');
		elseif SPMF
			[t,XYZ,QQ,U,k,s,w] = spm_projectionsF_ui('Results');
		end
		if isempty(t), spm_figure('Clear','Interactive'); return; end

		% get voxel sizes (in 'V') from mat file
		%-----------------------------------------------------------------------
		load([CWD '/SPM.mat']);
		dat{i} = struct('XYZ',XYZ,'t',t,'v',V(4:6));
	end

	spm_render(dat);
	return;
end



% Perform the rendering
%_______________________________________________________________________

set(spm_figure('FindWin','Interactive'),'Name','executing','Pointer','watch');

% Put the pre-computed data in the right format.
%---------------------------------------------------------------------------
load render
pre = cell(size(Matrixes,1),1);
for i=1:size(Matrixes,1)
	pre{i}=struct('MM',eval(Matrixes(i,:)),...
		'ren',eval(Rens(i,:)),...
		'dep',eval(Depths(i,:)),...
		'data',[],...
		'max',0);
	pre{i}.data = cell(size(dat,1),1);
end

spm_progress_bar('Init', size(dat,1)*size(pre,1), 'Making pictures', 'Number completed');

mx = zeros(size(pre,1),1)+eps;
for j=1:size(dat,1),

	% field should be sampled every mm.
	%---------------------------------------------------------------------------
	v =  dat{j}.v;
	t1       = [];
	XYZ1     = [];
	XYZ      = dat{j}.XYZ;
	xx       = XYZ(1,:);
	for x1=(-v(1)/2):(v(1)/2-eps*20)
		XYZ(1,:) = xx+x1;
		XYZ1 = [XYZ1 XYZ];
		t1   = [t1 dat{j}.t];
	end

	t0        = [];
	yy        = XYZ1(2,:);
	XYZ       = [];
	for y1=(-v(2)/2):(v(2)/2-eps*20)
		XYZ1(2,:) = yy+y1;
		XYZ       = [XYZ XYZ1];
		t0        = [t0 t1];
	end

	t        = [];
	XYZ1     = [];
	zz       = XYZ(3,:);
	for z1=(-v(3)/2):(v(3)/2-eps*20)
		XYZ(3,:)  = zz+z1;
		XYZ1      = [XYZ1 XYZ];
		t         = [t t0];
	end



	for i=1:size(pre,1)

		% transform from Taliarach space to space of the rendered image
		%---------------------------------------------------------------------------
		MM = pre{i}.MM;
		xyz      = (MM(1:3,1:3)*XYZ1);
		xyz(1,:) = xyz(1,:) + MM(1,4);
		xyz(2,:) = xyz(2,:) + MM(2,4);
		xyz(3,:) = xyz(3,:) + MM(3,4);
		xyz      = round(xyz);

		% only use values which will fit on the image
		%---------------------------------------------------------------------------
		dep = pre{i}.dep;

		msk = find((xyz(1,:) >= 1) & (xyz(1,:) <= size(dep,1)) ...
			&  (xyz(2,:) >= 1) & (xyz(2,:) <= size(dep,2)));

		X = zeros(size(dep));

		xyz      = xyz(:,msk);
		t0       = reshape(t(msk),1,prod(size(msk)));

		% calculate 'depth' of values
		%---------------------------------------------------------------------------
		z1  = dep(xyz(1,:)+(xyz(2,:)-1)*size(dep,1));
		msk = find(xyz(3,:) < (z1+60) & xyz(3,:) > (z1-5));

		if ~isempty(msk)

			% generate an image of the integral of the blob values.
			%---------------------------------------------------------------------------
			xyz = xyz(:,msk);
			dst = xyz(3,:) - z1(msk) - 5;
			dst = max(dst,0);
			t0  = t0(msk).*exp(-dst/10);
			X = full(sparse(xyz(1,:), xyz(2,:), t0, size(dep,1), size(dep,2)));
		end

		mx(j) = max([mx(j) max(max(X))]);
		pre{i}.data{j} = X;

		spm_progress_bar('Set', i+(j-1)*size(pre,1));
	end
end

mxmx = max(mx);
spm_progress_bar('Clear');
figure(spm_figure('FindWin','Graphics'));
spm_figure('Clear','Graphics');

% Combine the brain surface renderings with the blobs, and display using
% 24 bit colour.
%---------------------------------------------------------------------------
for i=1:size(Matrixes,1)
	ren = pre{i}.ren/64;
	X = cell(3,1);
	for j=1:size(pre{i}.data,1),
		X{j} = pre{i}.data{j}/mxmx;
	end
	for j=(size(pre{i}.data,1)+1):3
		X{j}=zeros(size(X{1}));
	end

	rgb = zeros([size(ren) 3]);
	tmp = ren.*max(1-X{1}-X{2}-X{3},0);
	rgb(:,:,1) = tmp + X{1};
	rgb(:,:,2) = tmp + X{2};
	rgb(:,:,3) = tmp + X{3};

	subplot(4,2,i);
	image(rgb);
	axis image off xy
end

spm_figure('Clear','Interactive');

