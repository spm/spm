function spm_render(dat,flg,rendfile)
% Render blobs on surface of a 'standard' brain
% FORMAT spm_render(dat,flg,rendfile)
%
% dat - a vertical cell array of length 1 to 3
%       - each element is a structure containing:
%         - XYZ - the x, y & z coordinates of the transformed t values.
%                 in units of voxels.
%         - t   - the SPM{.} values
%         - mat - affine matrix mapping from XYZ voxels to Talairach.
%         - dim - dimensions of volume from which XYZ is drawn.
% flg - optional.  If 1, then displays using the old style with hot
%       metal for the blobs, and grey for the brain.
% rendfile - the file containing the images to render on to. See also
%            spm_xbrain.m.
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

global SWD;

%=======================================================================
if nargin==0
	SPMid = spm('FnBanner',mfilename,'%I%');
	[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Results: render',0);

	rendfile = spm_get(1,'render*.mat',['Render file'],...
				fullfile(SWD,'rend'));

	num = spm_input('Number of sets',1,'1 set|2 sets|3 sets',[1 2 3]);
	flg = 0; if num==1, flg = spm_input('Style',1,'new|old',[0 1], 1); end;

	dat = cell(num,1);
	for i=1:num,
		[SPM,VOL] = spm_getSPM;
		% if isempty(SPM.Z),
		%	spm_figure('Clear',Finter)
		%	return
		% end
		dat{i} = struct(	'XYZ',	SPM.XYZ,...
					't',	SPM.Z',...
					'mat',	VOL.M,...
					'dim',	VOL.DIM);
	end

	spm_render(dat,flg,rendfile);
	return;
end

if nargin==1, flg = 0; end


% Perform the rendering
%=======================================================================
spm('Pointer','Watch')

load(rendfile);

if (exist('rend') ~= 1), % Assume old format...
	rend = cell(size(Matrixes,1),1);
	for i=1:size(Matrixes,1),
		rend{i}=struct('M',eval(Matrixes(i,:)),...
			'ren',eval(Rens(i,:)),...
			'dep',eval(Depths(i,:)));
		rend{i}.ren = rend{i}.ren/max(max(rend{i}.ren));
	end;
end;

spm_progress_bar('Init', size(dat,1)*length(rend),...
			'Formatting Renderings', 'Number completed');
for i=1:length(rend),
	rend{i}.max=0;
	rend{i}.data = cell(size(dat,1),1);
	if issparse(rend{i}.ren),
		% Assume that images have been DCT compressed
		% - the SPM99 distribution was originally too big.
		d = size(rend{i}.ren);
		B1 = spm_dctmtx(d(1),d(1));
		B2 = spm_dctmtx(d(2),d(2));
		rend{i}.ren = B1*rend{i}.ren*B2';
		% the depths did not compress so well with
		% a straight DCT - therefore it was modified slightly
		rend{i}.dep = exp(B1*rend{i}.dep*B2')-1;
	end;
	msk = find(rend{i}.ren>1);rend{i}.ren(msk)=1;
	msk = find(rend{i}.ren<0);rend{i}.ren(msk)=0;
	spm_progress_bar('Set', i);
end;
spm_progress_bar('Clear');

spm_progress_bar('Init', size(dat,1)*length(rend),...
			'Making pictures', 'Number completed');

mx = zeros(length(rend),1)+eps;
for j=1:size(dat,1),
	XYZ = dat{j}.XYZ;
	t   = dat{j}.t;
	dim = dat{j}.dim;
	mat = dat{j}.mat;

	for i=1:length(rend),

		% transform from Taliarach space to space of the rendered image
		%-------------------------------------------------------
		M1  = rend{i}.M*dat{j}.mat;
		zm  = sum(M1(1:2,1:3).^2,2).^(-1/2);
		M2  = diag([zm' 1 1]);
		M  = M2*M1;
		cor = [1 1 1 ; dim(1) 1 1 ; 1 dim(2) 1; dim(1) dim(2) 1 ;
		       1 1 dim(3) ; dim(1) 1 dim(3) ; 1 dim(2) dim(3); dim(1) dim(2) dim(3)]';
		tcor= M(1:3,1:3)*cor + M(1:3,4)*ones(1,8);
		off = min(tcor(1:2,:)');
		M2  = spm_matrix(-off+1)*M2;
		M  = M2*M1;
		xyz = (M(1:3,1:3)*XYZ + M(1:3,4)*ones(1,size(XYZ,2)));
		d2  = ceil(max(xyz(1:2,:)'));

		% calculate 'depth' of values
		%-------------------------------------------------------
		dep = spm_slice_vol(rend{i}.dep,spm_matrix([0 0 1])*inv(M2),d2,1);
		z1  = dep(round(xyz(1,:))+round(xyz(2,:)-1)*size(dep,1));

		if flg==1, msk = find(xyz(3,:) < (z1+20) & xyz(3,:) > (z1-5));
		else,      msk = find(xyz(3,:) < (z1+60) & xyz(3,:) > (z1-5)); end;

		if ~isempty(msk),

			% generate an image of the integral of the blob values.
			%-----------------------------------------------
			xyz = xyz(:,msk);
			if flg==1, t0  = t(msk);
			else,	dst = xyz(3,:) - z1(msk) - 5;
				dst = max(dst,0);
				t0  = t(msk).*exp(-dst/10)';
			end;
			X0  = full(sparse(xyz(1,:), xyz(2,:), t0, d2(1), d2(2)));
			hld = 1; if flg==1, hld = 0; end;
			X   = spm_slice_vol(X0,spm_matrix([0 0 1])*M2,size(rend{i}.dep),hld);
		else,
			X = zeros(size(rend{i}.dep));
		end;

		mx(j) = max([mx(j) max(max(X))]);
		rend{i}.data{j} = X;

		spm_progress_bar('Set', i+(j-1)*length(rend));
	end;
end;

mxmx = max(mx);

spm_progress_bar('Clear');
Fgraph = spm_figure('GetWin','Graphics');
spm_figure('Clear',Fgraph);

nrow = ceil(length(rend)/2);

if flg==1,
	% Old style split colourmap display.
	%---------------------------------------------------------------
	load Split;
	colormap(split);
	for i=1:length(rend),
		ren = rend{i}.ren;
		X   = rend{i}.data{1}/mxmx;
		msk = find(X);
		ren(msk) = X(msk)+1;
		subplot(nrow,2,i);
		image(ren*64);
		axis image off xy
	end;

else,
	% Combine the brain surface renderings with the blobs, and display using
	% 24 bit colour.
	%---------------------------------------------------------------
	for i=1:length(rend),
		ren = rend{i}.ren;
		X = cell(3,1);
		for j=1:size(rend{i}.data,1),
			X{j} = rend{i}.data{j}/mxmx;
		end
		for j=(size(rend{i}.data,1)+1):3
			X{j}=zeros(size(X{1}));
		end

		rgb = zeros([size(ren) 3]);
		tmp = ren.*max(1-X{1}-X{2}-X{3},0);
		rgb(:,:,1) = tmp + X{1};
		rgb(:,:,2) = tmp + X{2};
		rgb(:,:,3) = tmp + X{3};

		subplot(nrow,2,i);
		image(rgb);
		axis image off xy
	end;
end;

spm_figure('Clear','Interactive')
spm('Pointer');

