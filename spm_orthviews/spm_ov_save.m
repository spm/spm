function ret = spm_ov_save(varargin)
% Save as image tool - plugin for spm_orthviews
%
% This routine is a plugin to spm_orthviews. For general help about
% spm_orthviews and plugins type
%             help spm_orthviews
% at the MATLAB prompt.
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_ov_save.m 5534 2013-06-10 13:29:09Z guillaume $


cmd = lower(varargin{1});
switch cmd
    % Context menu and callbacks
    case 'context_menu'
        item0 = uimenu(varargin{3}, ...
            'Label', 'Save As...', ...
            'Callback', @orthviews_save);
        ret = item0;
    otherwise
end

%==========================================================================
function orthviews_save(hObj,event)

global st

hC = current_handle;

[filename, pathname, filterindex] = uiputfile(...
    {'*.png' 'PNG files (*.png)'}, 'Save as');
if isequal(filename,0) || isequal(pathname,0), return; end

imgfile = fullfile(pathname, filename);

X  = frame2im(getframe(st.fig));
sz = size(X);
sz = [sz(1) sz(1) sz(2) sz(2)];
p1 = get(st.vols{hC}.ax{1}.ax,'Position');
p2 = get(st.vols{hC}.ax{3}.ax,'Position');
a  = [p1(1) p1(2)  p2(1)+p2(3)-p1(1) p2(2)+p2(4)-p1(2)] + 0.005*[-1 -1 2 2];
a  = max(min(a,1),0);
sz = ([1-a(2)-a(4),1-a(2),a(1),a(1)+a(3)] .* (sz-1)) + 1;
sz = round(sz);
X  = X(sz(1):sz(2),sz(3):sz(4),:);
imwrite(X,imgfile,'png');

%==========================================================================
function h = current_handle
global st
try
    hs = [];
    for i = spm_orthviews('valid_handles')
        hs = [hs st.vols{i}.ax{1}.cm];
    end
    hc = get(gca,'UIContextMenu');
    h  = find(hs==hc);
catch
    h  = 1;
end
