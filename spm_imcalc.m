function sf=spm_imcalc(Vi,Vo,f,hold,mask)
% Perform algebraic functions on images
% FORMAT sf=spm_imcalc(Vi,Vo,f,hold,mask)
% Vi    - vector of mapped image volumes to work on (from spm_vol)
% Vo    - handle structure for output image volume
% f     - expression to be evaluated
% hold  - interpolation hold (see spm_slice_vol)
%         [defaults (missing or empty) to 0 - nearest neighbour]
% mask  - implicit zero mask?
%         [defaults (missing or empty) to 0]
% sf    - scalefactor for output image
%
% With no arguments, spm_imcalc_ui is called.
%_______________________________________________________________________
%
% The images specified in Vi, are referred to as i1, i2, i3,...
% in the expression to be evaluated.
%
% See spm_imcalc_ui for example usage...
%_______________________________________________________________________
% %W% John Ashburner, Andrew Holmes %E%


%-Parameters & arguments
%=======================================================================
if nargin<5, mask=[]; end, if isempty(mask), mask=0; end
if nargin<4, hold=[]; end, if isempty(hold), hold=0; end
if nargin==0, sf=spm_imcalc_ui; return, end
if nargin<3, error('insufficient arguments'), end

if Vo.dim(4)~=4, error('only 16 bit signed short output images supported'), end


%=======================================================================
%-Computation
%=======================================================================
n   = prod(size(Vi));			%-#images
if n==0, error('no input images specified'), end
Y   = zeros(Vo.dim(1:3));		%-result of calculations


%-Start progress plot
%-----------------------------------------------------------------------
spm_progress_bar('Init',Vo.dim(3),f,'planes completed');


%-Loop over planes computing result Y
%-----------------------------------------------------------------------
for p = 1:Vo.dim(3),
	B = spm_matrix([0 0 -p 0 0 0 1 1 1]);

	for i = 1:n
		M = inv(B*inv(Vo.mat)*Vi(i).mat);
		d = spm_slice_vol(Vi(i),M,Vo.dim(1:2),hold);
		if mask & Vi.dim(4)<16, d(d==0)=NaN; end
		eval(['i',num2str(i) '=d;']);
	end

	eval(['yp = ' f ';'],['error([''Can''t evaluate "'',f,''".'']);']);
	if (prod(Vo.dim(1:2)) ~= prod(size(yp)))
		error(['"',f,'" produced incompatible image.']); end
	Y(:,:,p) = yp;

	spm_progress_bar('Set',p);
end


%-Write output image (16 bit signed)
%=======================================================================
sf = max(abs(Y(:)))/32767;		%-compute scale factor
if sf==0, sf = 1; end			%-for images of all zeros
Vo.pinfo = [sf 0 0]';			%-set scalefactor
Y(isnan(Y))=0;				%-implicit zero mask for spm_type(4)
for p = 1:Vo.dim(3), spm_write_plane(Vo,Y(:,:,p),p); end


%-End
%-----------------------------------------------------------------------
spm_progress_bar('Clear')
