function spm_check_registration(varargin)
% A visual check of image registration quality
% FORMAT spm_check_registration
% FORMAT spm_check_registration(images)
% Orthogonal views of one or more images are displayed. Clicking in
% any image moves the centre of the orthogonal views. Images are
% shown in orientations relative to that of the first selected image.
% The first specified image is shown at the top-left, and the last at
% the bottom right. The fastest increment is in the left-to-right
% direction (the same as you are reading this).
%__________________________________________________________________________
% Copyright (C) 1997-2013 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_check_registration.m 5700 2013-10-17 14:59:50Z guillaume $

SVNid = '$Rev: 5700 $';

%-Get input
%--------------------------------------------------------------------------
if nargin
    if nargin > 1 && iscellstr(varargin{2})
        images = varargin{1};
    elseif isstruct(varargin{1})
        images = [varargin{:}];
    else
        images = char(varargin);
    end
else
    [images, sts] = spm_select([1 24],'image','Select images');
    if ~sts, return; end
end

if ischar(images), images = spm_vol(images); end
images = images(1:min(numel(images),24));

%-Print
%--------------------------------------------------------------------------
spm('FnBanner',mfilename,SVNid);                                        %-#
exactfname  = @(f) [f.fname ',' num2str(f.n(1))];
if desktop('-inuse')
    href    = '<a href="matlab:%s;">%s</a>';
    cmd     = 'spm_image(''display'',''%s'')';
    dispone = @(f) sprintf(href,sprintf(cmd,f),f);
    cmd     = 'spm_check_registration(%s)';
    str     = [];
    for i=1:numel(images)
        str = [str sprintf('''%s'',',exactfname(images(i)))];
    end
    dispall = @(f) sprintf([' (' href ')  '],sprintf(cmd,str(1:end-1)),'all');
else
    dispone = @(f) f;
    dispall = @(f) '        ';
end
for i=1:numel(images)
    if i==1,     fprintf('Display ');                                   %-#
    elseif i==2, fprintf('%s',dispall(images));
    else         fprintf('        '); end
    fprintf('%s\n',dispone(exactfname(images(i))));                     %-#
end

%-Display
%--------------------------------------------------------------------------
spm_figure('GetWin','Graphics');
spm_figure('Clear','Graphics');
spm_orthviews('Reset');

mn = length(images);
n  = round(mn^0.4);
m  = ceil(mn/n);
w  = 1/n;
h  = 1/m;
ds = (w+h)*0.02;
for ij=1:mn
    i = 1-h*(floor((ij-1)/n)+1);
    j = w*rem(ij-1,n);
    handle = spm_orthviews('Image', images(ij),...
        [j+ds/2 i+ds/2 w-ds h-ds]);
    if ij==1, spm_orthviews('Space'); end
    spm_orthviews('AddContext',handle);
end

%-Backward compatibility with spm_check_registration(images,captions)
%--------------------------------------------------------------------------
if nargin > 1 && iscellstr(varargin{2})
    spm_orthviews('Caption',varargin{2}, varargin{3:end});
    %for ij=1:numel(varargin{2})
    %    spm_orthviews('Caption', ij, varargin{2}{ij}, varargin{3:end});
    %end
end
