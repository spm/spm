function spm_mip(X,L,V)
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
% %W% Karl Friston %E%


%-Get GRID value
%-----------------------------------------------------------------------
global GRID, if isempty(GRID), GRID = 0.6; end

% single slice case (ORIGIN = [0 0])
%=======================================================================
if V(3) == 1
	i   = L(1,:)/V(4);
	j   = L(2,:)/V(5);
	d   = i > 1 & i < V(1) & j > 1 & j < V(2);
	mip = full(sparse(i(d),j(d),X(d),V(1),V(2)));
	imagesc([1 V(1)*V(4)],[1 V(2)*V(5)],(-mip')); axis xy
	axis image; 
	set(gca,'FontSize',8,'TickDir','in')
	xlabel('x'), ylabel('y')
	return
end

%-3d case
%=======================================================================

%-Remove negtive values from point list and scale to a maximium of unity
%-----------------------------------------------------------------------
X    = X(:)';
d    = X > 0;
L    = round(L(:,d));
X    = X(d);
X    = X/max(X);

%-Load mip and create maximum intensity projection
%-----------------------------------------------------------------------
load MIP
mip  = mip96*GRID;
d    = zeros(size(mip));
spm_project(X,L,d,V(1:6));
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
