function opts = spm_config_checkreg
% Configuration file for check-reg jobs
%_______________________________________________________________________
% %W% %E%

%_______________________________________________________________________

w = spm_jobman('HelpWidth');

%_______________________________________________________________________

data.type = 'files';
data.name = 'Images to Display';
data.tag  = 'data';
data.filter = 'image';
data.num  = Inf;
data.help = {'Images to display.'};

opts.type = 'branch';
opts.name = 'Check Registration';
opts.tag  = 'checkreg';
opts.val  = {data};
opts.prog = @dispims;
opts.help = spm_justify(w,[...
' Orthogonal views of one or more images are displayed.  Clicking in'...
' any image moves the centre of the orthogonal views.  Images are'...
' shown in orientations relative to that of the first selected image.'...
' The first specified image is shown at the top-left, and the last at'...
' the bottom right.  The fastest increment is in the left-to-right'...
' direction (the same as you are reading this).']);
%------------------------------------------------------------------------

%------------------------------------------------------------------------
function dispims(varargin)
job = strvcat(varargin{1}.data);
spm_check_registration(job);
return;
