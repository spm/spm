function ret = spm_ov_display(varargin)
% Display tool - plugin for spm_orthviews
%
% This routine is a plugin to spm_orthviews. For general help about
% spm_orthviews and plugins type
%             help spm_orthviews
% at the MATLAB prompt.
%__________________________________________________________________________
% Copyright (C) 2013-2014 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_ov_display.m 6340 2015-02-16 12:25:56Z guillaume $


switch lower(varargin{1})
    % Context menu and callbacks
    case 'context_menu'
        item0 = uimenu(varargin{3}, ...
            'Label', 'Display');
        item1 = uimenu(item0, ...
            'Label', 'Intensities', ...
            'Tag', 'OVmenu_Intensities',...
            'Callback', @orthviews_display);
        item2 = uimenu(item0, ...
            'Label', 'Filenames', ...
            'Tag', 'OVmenu_Filenames',...
            'Callback', @orthviews_display);
        item3 = uimenu(item0, ...
            'Label', 'Coordinates', ...
            'Tag', 'OVmenu_Coordinates',...
            'Callback', @orthviews_display);
        item4 = uimenu(item0, ...
            'Label', 'Labels');
        list = spm_atlas('List','installed');
        for i=1:numel(list)
            uimenu(item4, ...
            'Label', list(i).name, ...
            'Tag', ['OVmenu_' list(i).name],...
            'Callback', @orthviews_display);
        end
        ret = item0;
    case 'redraw'
        orthviews_display_redraw(varargin{2:end});
    case 'label'
        ret = 'Display';
    otherwise
end

%==========================================================================
function orthviews_display(hObj,event)

global st

if strcmp(get(hObj, 'Checked'),'on')
    set(findobj(st.fig,'-regexp','Tag','^OVmenu_'), 'Checked', 'off');
    for i=1:numel(st.vols)
        if ~isempty(st.vols{i})
            xlabel(st.vols{i}.ax{3}.ax,'');
            try, st.vols{i} = rmfield(st.vols{i},'display'); end
        end
    end
else
    set(findobj(st.fig,'-regexp','Tag','^OVmenu_'), 'Checked', 'off');
    set(findobj(st.fig,'Tag',['OVmenu_' get(hObj,'Label')]), 'Checked', 'on');
    dsp = get(hObj,'Label');
    if ~ismember(dsp,{'Intensities','Filenames','Coordinates'})
        dsp = spm_atlas('Load',dsp);
    end
    for i=1:numel(st.vols)
        if ~isempty(st.vols{i})
            st.vols{i}.display = dsp;
        end
    end
    spm_ov_display('redraw');
end

%==========================================================================
function orthviews_display_redraw(i,varargin) %i, TM0, TD, CM0, CD, SM0, SD

global st

if ~nargin, n = 1:numel(st.vols); else n = i; end

for i=n
    if ~isempty(st.vols{i})
        action = st.vols{i}.display;
        if ~ischar(action), xA = action; action = 'Labels'; end
        switch action
            case 'Intensities'
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
            case 'Filenames'
                Ys = [spm_file(st.vols{i}.fname,'filename') ',' num2str(st.vols{i}.n(1))];
            case 'Coordinates'
                XYZ   = spm_orthviews('pos',i);
                XYZmm = st.vols{i}.mat(1:3,:)*[XYZ;1];
                Ys    = sprintf('mm: %0.1f %0.1f %0.1f\nvx: %0.1f %0.1f %0.1f',XYZmm,XYZ);
            case 'Labels'
                pos = spm_orthviews('pos',i);
                try
                    Ys  = spm_atlas('Query',xA,st.vols{i}.mat(1:3,:)*[pos;1]);
                    Ys  = strrep(Ys,'_',' ');
                    Ys  = sprintf(strrep(Ys,'.','\n'));
                catch
                    Ys  = 'Labels';
                end
            otherwise
                Ys = 'error';
        end
        xlabel(st.vols{i}.ax{3}.ax,Ys,'Interpreter','none');
    end
end
