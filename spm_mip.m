function spm_mip(X,L,V)
% SPM maximum intensity projection
% FORMAT spm_mip(X,L,V);
% V  -  SPM
% L  -  Talairach coordinates
% V  -  {1 x 6} vector of image and voxel sizes [DIM VOX]
%___________________________________________________________________________
%
% If the data are 2 dimensional [V(3) = 1] the projection is simply an
% image, otherwise:
%
% spm_mip creates and displays a maximum intensity projection of a point
% list of voxel values (X) and their location (L) in three orthogonal
% views of the brain.  It is assumed voxel locations conform to the space
% defined in the atlas of Talairach and Tournoux (1988), otherwise disable
% the GRID.
%
% This routine loads mip (in MIP.mat).  mip is a 360 x 352 matrix with
% contours and grids defining the space of Talairach & Tournoux (1988)
% A default colormap of 64 levels is assumed.
%
% If global XVIEWER is set to a PGM viewer program, images are also
% displayed with that program.
%
%__________________________________________________________________________
% %W% %E%

global GRID XVIEWER TWD

% default GRID value
%---------------------------------------------------------------------------
if isempty(GRID)
	GRID = 0.6; end


% single slice case [ORIGIN = (0.0)]
%---------------------------------------------------------------------------
if V(3) == 1
	i   = L(1,:)/V(4);
	j   = L(2,:)/V(5);
	d   = i > 1 & i < V(1) & j > 1 & j < V(2);
	mip = full(sparse(i(d),j(d),X(d),V(1),V(2)));
	imagesc([1 V(1)*V(4)],[1 V(2)*V(5)],(mip')); axis xy
	axis off; axis image; 
	return
end

% remove negtive values from point list and scale to a maximium of unity
%---------------------------------------------------------------------------
X    = X(:)';
d    = X > 0;
L    = round(L(:,d));
X    = X(d);
X    = X/max(X);

% load mip and create maximum intensity projection
%---------------------------------------------------------------------------
load MIP
mip  = mip*GRID;
d    = zeros(size(mip));
spm_project(X,L,d,V);
mip  = max(d,mip);
image(rot90((1 - mip)*64)); axis image; axis off;

%
%
%
if ~isempty(XVIEWER) & isstr(XVIEWER)
        if isempty(TWD); TWD = '/tmp'; end
        % This should really be a subroutine...
        t = clock;
        TmpNm = sprintf('/%s/mip%02d%02d%02d%02d%02d.pgm',TWD, ...
                                  floor(t(3:6)),floor(100*(t(6)-floor(t(6)))));
        fid = fopen(TmpNm,'w');
        fprintf(fid,'P5\n%d %d\n255\n',size(mip));
        fwrite(fid,fliplr(mip*255),'uchar');
        fclose(fid);
        unix(['(' XVIEWER ' ' TmpNm ';\rm ' TmpNm ' ) &']);
end
