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
% Copyright (C) 1997-2014 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_check_registration.m 6245 2014-10-15 11:22:15Z guillaume $

SVNid = '$Rev: 6245 $';

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
if numel(images) > 24
    if ~isdeployed, addpath(fullfile(spm('Dir'),'spm_orthviews')); end
    img = cell(1,numel(images));
    for i=1:numel(images)
        img{i} = [images(i).fname ',' num2str(images(i).n(1))];
    end
    spm_ov_browser('ui',char(img));
    return
end
images = images(1:min(numel(images),24));

%-Print
%--------------------------------------------------------------------------
spm('FnBanner',mfilename,SVNid);                                        %-#
exactfname  = @(f) [f.fname ',' num2str(f.n(1))];
cmddispone  = 'spm_image(''display'',''%s'')';
cmddispall  = 'spm_check_registration(''%s'')';
if spm_platform('desktop')
    str     = '';
    for i=1:numel(images)
        str = [str sprintf('''%s'',',exactfname(images(i)))];
    end
    dispall = [' (' spm_file('all','link',sprintf(cmddispall,str(2:end-2))) ')  '];
else
    dispall = '        ';
end
for i=1:numel(images)
    if i==1,     fprintf('Display ');                                   %-#
    elseif i==2, fprintf('%s',dispall);
    else         fprintf('        '); end
    fprintf('%s\n',spm_file(exactfname(images(i)),'link',cmddispone));  %-#
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
