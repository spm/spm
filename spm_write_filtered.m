function Vo = spm_write_filtered(Z,XYZ,M,DIM,descrip)
% Writes the filtered SPM as an image
% FORMAT spm_write_filtered(Z,XYZ,M,DIM,descrip)
%
% Z       - {1 x ?} vector point list of SPM values for MIP
% XYZ     - {3 x ?} matrix of coordinates of points (Talairach coordinates)
% M       - voxels - > mm matrix
% DIM     - image dimensions {voxels}
% descrip - description string [default 'SPM-filtered']
%
%-----------------------------------------------------------------------
%
% spm_write_filtered takes a pointlist image (parallel matrixes of
% co-ordinates and voxel intensities), and writes it out into an image
% file.
%
% It is intended for writing out filtered SPM's from the results
% section of SPM, but can be used freestanding.
%_______________________________________________________________________
% %W% FIL %E%

%-Parse arguments
%-----------------------------------------------------------------------
if nargin<5, descrip='SPM-filtered'; end
if nargin<4, error('Insufficient arguments'), end

%-Get filename
%-----------------------------------------------------------------------
Q       = spm_str_manip(spm_input('Output filename',1,'s'),'sdv');
spm('Pointer','Watch')

%-Set up header information
%-----------------------------------------------------------------------
Vo      = struct(...
		'fname',	Q,...
		'dim',		[DIM', spm_type('uint8')],...
		'mat',		M,...
		'descrip', 	descrip);

%-Reconstruct filtered image from XYZ & Z
%-----------------------------------------------------------------------
Y      = zeros(DIM(1:3)');
iM     = inv(M);
rcp    = round(iM(1:3,:)*[XYZ; ones(1,size(XYZ,2))]);
OFF    = rcp(1,:) + DIM(1)*(rcp(2,:) + DIM(2)*rcp(3,:));
Y(OFF) = Z.*(Z > 0);

%-Write the filtered volume
%-----------------------------------------------------------------------
Vo = spm_write_vol(Vo,Y);
fprintf('\n%s: %s\n\n',mfilename,spm_get('CPath',Q,pwd))

%-End
%-----------------------------------------------------------------------
spm('Pointer','Arrow');
