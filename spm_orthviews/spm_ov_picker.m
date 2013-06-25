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
% $Id: spm_ov_picker.m 5564 2013-06-25 15:20:27Z guillaume $


switch lower(varargin{1})
    % Context menu and callbacks
    case 'context_menu'
        item0 = uimenu(varargin{3}, ...
            'Label', spm_ov_picker('label'), ...
            'Callback', @orthviews_picker);
        ret = item0;
    case 'redraw'
        orthviews_picker_redraw(varargin{2:end});
    case 'label'
        ret = 'Display intensities';
    otherwise
end

%==========================================================================
function orthviews_picker(hObj,event)

global st

if strcmp(get(hObj, 'Checked'),'on')
    set(findobj(st.fig,'Type','uimenu','Label',spm_ov_picker('label')), ...
        'Checked', 'off');
    for i=1:numel(st.vols)
        if ~isempty(st.vols{i})
            xlabel(st.vols{i}.ax{3}.ax,'');
            st.vols{i} = rmfield(st.vols{i},'picker');
        end
    end
else 
    set(findobj(st.fig,'Type','uimenu','Label',spm_ov_picker('label')), ...
        'Checked', 'on');
    for i=1:numel(st.vols)
        if ~isempty(st.vols{i})
            st.vols{i}.picker = [];
        end
    end
    spm_ov_picker('redraw');
end

%==========================================================================
function orthviews_picker_redraw(i,varargin) %i, TM0, TD, CM0, CD, SM0, SD

global st

if ~nargin, n = 1:numel(st.vols); else n = i; end

for i=n
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
