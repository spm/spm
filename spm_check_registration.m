function spm_check_registration(images)
% A visual check of image registration quality.
% FORMAT spm_check_registration
% Orthogonal views of two images are displayed.  Clicking in any image
% moves the centre of the orthogonal views.
%_______________________________________________________________________
% %W% John Ashburner %E%

if nargin==0
	images=spm_get(2,'.img',['Select two images']);
	spm_check_registration(images);
elseif nargin==1
	fg = spm_figure('Findwin','Graphics');
	if isempty(fg)
		fg=spm_figure('Create','Graphics');
		if isempty(fg)
			error('Cant create graphics window');
		end
	else
		spm_figure('Clear','Graphics');
	end
	h1=spm_orthviews('Image', images(1,:), [0.0 0.0 1 0.5]);
	spm_orthviews('Space',1);
	spm_orthviews('Image', images(2,:), [0.0 0.5 1 0.5]);
else
	error('Incorrect Usage');
end
