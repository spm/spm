function spm_render(dat)
% Render blobs on surface of a 'standard' brain
% FORMAT spm_render
%_______________________________________________________________________
% 
% FORMAT spm_render(dat)
% dat - a vertical cell array of length 1 to 3
%       - each element is a structure containing:
%         - XYZ - the x, y & z coordinated of the transformed t values.
%                 in units of voxels.
%         - t   - the transformed t values
%         - mat - affine matrix mapping from XYZ to Talairach.
%         - dim - dimensions of volume from which XYZ is drawn.
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


if nargin==0,

	spm_figure('Clear','Interactive');
	Finter = spm_figure('FindWin','Interactive');
	set(Finter,'Name','SPM render')

	num = spm_input('Number of sets',2,'1 set|2 sets|3 sets',[1 2 3]);
	dat = cell(num,1);
	for i=1:num,
		[SPM,VOL,DES] = spm_getSPM;
		% if isempty(SPM.Z), spm_figure('Clear','Interactive'); return; end
		imat   = inv(VOL.M);
		XYZ    = round(imat(1:3,1:3)*VOL.XYZ + imat(1:3,4)*ones(1,size(VOL.XYZ,2)));
		dat{i} = struct('XYZ',XYZ, 't',SPM.Z', 'mat',VOL.M, 'dim', VOL.DIM);
	end;

	spm_render(dat);
	return;
end;

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
	XYZ = dat{j}.XYZ;
	t   = dat{j}.t;
	dim = dat{j}.dim;
	mat = dat{j}.mat;

	for i=1:size(pre,1)

		% transform from Taliarach space to space of the rendered image
		%---------------------------------------------------------------------------
		M1  = pre{i}.MM*dat{j}.mat;
		zm  = sum(M1(1:2,1:3).^2,2).^(-1/2);
		M2  = diag([zm' 1 1]);
		MM  = M2*M1;
		cor = [1 1 1 ; dim(1) 1 1 ; 1 dim(2) 1; dim(1) dim(2) 1 ;
		       1 1 dim(3) ; dim(1) 1 dim(3) ; 1 dim(2) dim(3); dim(1) dim(2) dim(3)]';
		tcor= MM(1:3,1:3)*cor + MM(1:3,4)*ones(1,8);
		off = min(tcor(1:2,:)');
		M2  = spm_matrix(-off+1)*M2;
		MM  = M2*M1;
		xyz = (MM(1:3,1:3)*XYZ + MM(1:3,4)*ones(1,size(XYZ,2)));
		d2  = ceil(max(xyz(1:2,:)'));

		% calculate 'depth' of values
		%---------------------------------------------------------------------------
		dep = spm_slice_vol(pre{i}.dep,spm_matrix([0 0 1])*inv(M2),d2,1);
		z1  = dep(round(xyz(1,:))+round(xyz(2,:)-1)*size(dep,1));
		msk = find(xyz(3,:) < (z1+60) & xyz(3,:) > (z1-5));

		if ~isempty(msk),

			% generate an image of the integral of the blob values.
			%---------------------------------------------------------------------------
			xyz = xyz(:,msk);
			dst = xyz(3,:) - z1(msk) - 5;
			dst = max(dst,0);
			t0  = t(msk).*exp(-dst/10)';
			X0  = full(sparse(xyz(1,:), xyz(2,:), t0, d2(1), d2(2)));
			X   = spm_slice_vol(X0,spm_matrix([0 0 1])*M2,size(pre{i}.dep),1);
		else,
			X = zeros(size(dep));
		end;

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

