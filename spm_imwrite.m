function spm_imwrite(V,dat)
% Write an image volume to disk.
% FORMAT spm_imwrite(V,dat)
% V   - a structure containing image volume information (see spm_vol).
% dat - a one, two or three dimensional matrix containing the image voxels.
% ____________________________________________________________________________
% %W% John Ashburner %E%

dim = size(dat);
if length(dim)>3,
	error('Can only handle a maximum of 3 dimensions.');
end;
dim = [dim 1 1 1];
if ~all(dim(1:3) == V.dim(1:3)) | (size(V.pinfo,2)~=1 & size(V.pinfo,2)~=dim(3)),
	error('Incompatible dimensions.');
end;


% Set scalefactors and offsets
dt = V.dim(4); if dt>256, dt = dt/256; end;
s = find(dt == [2 4 8 128+2 128+4 128+8]);
datminmax = [0 -2^15 -2^31 0 0 ; 2^8-1 2^15-1 2^31-1 2^16 2^32];
datminmax = datminmax(:,s);
V.pinfo(1,:)=1;
V.pinfo(2,:)=0;
if ~isempty(s),
	if size(V.pinfo,2) ~= 1,
		for p=1:dim(3),
			mx = max(max(dat(:,:,p)));
			mn = min(min(dat(:,:,p)));
			V.pinfo(1,p) = (mx-mn)/(datminmax(2)-datminmax(1));
			V.pinfo(2,p) = mx-datminmax(1)/V.pinfo(1,p);
		end;
	else,
		mx = max(max(max(dat)));
		mn = min(min(min(dat)));
		V.pinfo(1,1) = (mx-mn)/(datminmax(2)-datminmax(1));
		V.pinfo(2,1) = mx-datminmax(1)/V.pinfo(1,1);
	end;
end;

spm_create_image(V);
for p=1:V.dim(3),
	spm_write_plane(V,dat(:,:,p),p);
end;

return;
