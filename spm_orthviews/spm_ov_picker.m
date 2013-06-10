function ret = spm_ov_picker(varargin)
% Picker tool - plugin for spm_orthviews
%
% This routine is a plugin to spm_orthviews. For general help about
% spm_orthviews and plugins type
%             help spm_orthviews
% at the MATLAB prompt.
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_ov_picker.m 5534 2013-06-10 13:29:09Z guillaume $


cmd = lower(varargin{1});
switch cmd
    % Context menu and callbacks
    case 'context_menu'
        item0 = uimenu(varargin{3}, ...
            'Label', 'Display values', ...
            'Callback', @orthviews_picker);
        ret = item0;
    case 'update'
        orthviews_picker_update;
    otherwise
end

%==========================================================================
function orthviews_picker(hObj,event)

global st

if strcmp(get(hObj, 'Checked'),'on')
    set(hObj, 'Checked', 'off');
    st.callback = ';';
    for i=1:numel(st.vols)
        if ~isempty(st.vols{i})
            xlabel(st.vols{i}.ax{3}.ax,'');
        end
    end
else 
    set(hObj, 'Checked', 'on');
    st.callback = 'spm_ov_picker(''update'');';
    eval(st.callback);
end

%==========================================================================
function orthviews_picker_update

global st

for i=1:numel(st.vols)
    if ~isempty(st.vols{i})
        pos = spm_orthviews('pos',i);
        try
            Y = spm_sample_vol(st.vols{i},pos(1),pos(2),pos(3),st.hld);
        catch
            Y = NaN;
            fprintf('Cannot access file "%s".\n', st.vols{i}.fname);
        end
        Ys = sprintf('Y = %g',Y);
        %Y = get(findobj(st.vols{i}.ax{1}.cm,'UserData','v_value'),'Label');
        if isfield(st.vols{i},'blobs')
            for j=1:numel(st.vols{i}.blobs)
                p = st.vols{i}.blobs{j}.mat\st.vols{i}.mat*[pos;1];
                Y = spm_sample_vol(st.vols{i}.blobs{j}.vol,p(1),p(2),p(3),st.hld);
                Ys = sprintf('%s\nY = %g',Ys,Y);
            end
        end
        xlabel(st.vols{i}.ax{3}.ax,Ys);
    end
end
