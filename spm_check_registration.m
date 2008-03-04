function spm_check_registration(images)
% A visual check of image registration quality.
% FORMAT spm_check_registration
% Orthogonal views of one or more images are displayed.  Clicking in
% any image moves the centre of the orthogonal views.  Images are
% shown in orientations relative to that of the first selected image.
% The first specified image is shown at the top-left, and the last at
% the bottom right.  The fastest increment is in the left-to-right
% direction (the same as you are reading this).
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_check_registration.m 1185 2008-03-04 16:31:21Z volkmar $


if nargin==0,
    images = spm_select([1 15],'image','Select images');
    if size(images,1)<1, return; end;
    spm_check_registration(images);
elseif nargin==1,
    fg = spm_figure('Findwin','Graphics');
    if isempty(fg),
        fg=spm_figure('Create','Graphics');
        if isempty(fg),
            error('Cant create graphics window');
        end;
    else
        spm_figure('Clear','Graphics');
    end;
    if isstruct(images) && isfield(images,'data')
        images = spm_vol(strvcat(images.data));
    elseif ischar(images), images=spm_vol(images); end;
    spm_orthviews('Reset');
    mn = length(images);
    n  = round(mn^0.4);
    m  = ceil(mn/n);
    w  = 1/n;
    h  = 1/m;
    ds = (w+h)*0.02;
    for ij=1:mn,
        i  = 1-h*(floor((ij-1)/n)+1);
        j  = w*rem(ij-1,n);
        handle(ij) = spm_orthviews('Image', images(ij),...
            [j+ds/2 i+ds/2 w-ds h-ds]);
        if ij==1, spm_orthviews('Space'); end;
        spm_orthviews('AddContext',handle(ij));
    end;
else
    error('Incorrect Usage');
end;
