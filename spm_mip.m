function spm_mip(X,VOL)
% SPM maximum intensity projection
% FORMAT spm_mip(X,L,V);
% V  -  SPM
% L  -  Talairach coordinates
% V  -  {6 x 1} vector of image & voxel sizes [DIM VOX]
%    or {9 x 1} vector of image & voxel sizes, and origin [DIM VOX ORIGIN]
%       (ORIGIN is required for 2D MIPs)
%_______________________________________________________________________
%
% If the data are 2 dimensional [V(3) = 1] the projection is simply an
% image, otherwise:
%
% spm_mip creates and displays a maximum intensity projection of a point
% list of voxel values (X) and their location (L) in three orthogonal
% views of the brain.  It is assumed voxel locations conform to the space
% defined in the atlas of Talairach and Tournoux (1988).
%
% This routine loads a mip putline from MIP.mat. This is an image with
% contours and grids defining the space of Talairach & Tournoux (1988).
% mip95 corresponds to the Talairach atlas, mip96 to the MNI templates.
% The outline and grid are superimposed at intensity GRID (global default),
% defaulting to 0.6.
%
% A default colormap of 64 levels is assumed.
%
% If global XVIEWER is set to a PGM viewer program, images are also
% displayed with that program.
%
%_______________________________________________________________________
% %W% Karl Friston and others %E%


%-Get GRID value
%-----------------------------------------------------------------------
global GRID, if isempty(GRID), GRID = 0.6; end

%-Remove negative values from point list and scale to a maximium of unity
%-----------------------------------------------------------------------
X    = X(:)';
d    = find(X > 0);
L    = VOL.XYZ(:,d);
X    = X(d);
X    = X/max(X);

% single slice case (ORIGIN = [0 0])
%=======================================================================
if VOL.DIM(3) == 1,
	vox = sqrt(sum(VOL.M(1:3,1:3).^2));
	L   = round(VOL.M\[L ; ones(1,size(L,2))]);
	mip = full(sparse(L(1,:),L(2,:),X,VOL.DIM(1),VOL.DIM(2)));
	imagesc([1 VOL.DIM(1)*vox(1)],[1 VOL.DIM(2)*vox(2)],-mip');
	axis xy image; 
	set(gca,'FontSize',8,'TickDir','in')
	xlabel('x'); ylabel('y');
	return;
end;

%-3d case
%=======================================================================
%-Load mip and create maximum intensity projection
%-----------------------------------------------------------------------
load MIP
mip  = mip96*GRID;
c    = [0 0 0 ; 0 0 1 ; 0 1 0 ; 0 1 1 
	1 0 0 ; 1 0 1 ; 1 1 0 ; 1 1 1]-0.5;
c    = (VOL.M(1:3,1:3)*c')';
dim  = [(max(c)-min(c)) size(mip)];
d    = spm_project(X,round(L),dim);
mip  = max(d,mip);
image(rot90((1 - mip)*64)); axis image; axis off;

%-PGM file to viewer app
%=======================================================================
global XVIEWER TWD
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
