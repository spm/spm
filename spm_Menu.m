function Fmenu = spm_Menu(action, varargin)
% SPM Menu window
% FORMAT Fmenu = spm_Menu('Create',Vis)
% Create the SPM Menu window (tag property set to 'Menu')
% Vis    - visibility: {'on'} or 'off'
% Fmenu  - handle of figure created
%
% FORMAT Fmenu = spm_Menu('Switch',Modality)
% Switch the SPM Menu window to the specified modality
%
% FORMAT spm_Menu('Close')
% Close the SPM Menu window
%_________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_Menu.m 7432 2018-09-28 14:35:27Z guillaume $


if nargin < 1, action = 'Create'; end

switch lower(action)
    case 'create'
        spm_Menu('Close');
        Fmenu = create_menu(varargin{1:end});
        
    case 'switch'
        Fmenu = spm_figure('FindWin','Menu');
        [Modality,ModNum] = spm('CheckModality',varargin{1});
        if ~isempty(Fmenu)
            if strcmpi(Modality,'PET')
                set(findobj(Fmenu, '-regexp', 'Tag', 'FMRI'), 'Visible', 'off');
                set(findobj(Fmenu, '-regexp', 'Tag', 'EEG'),  'Visible', 'off');
                set(findobj(Fmenu, '-regexp', 'Tag', 'PET'),  'Visible', 'on' );
            elseif strcmpi(Modality,'FMRI')
                set(findobj(Fmenu, '-regexp', 'Tag', 'EEG'),  'Visible', 'off');
                set(findobj(Fmenu, '-regexp', 'Tag', 'PET'),  'Visible', 'off');
                set(findobj(Fmenu, '-regexp', 'Tag', 'FMRI'), 'Visible', 'on' );
            else
                set(findobj(Fmenu, '-regexp', 'Tag', 'PET'),  'Visible', 'off');
                set(findobj(Fmenu, '-regexp', 'Tag', 'FMRI'), 'Visible', 'off');
                set(findobj(Fmenu, '-regexp', 'Tag', 'EEG'),  'Visible', 'on' );
            end
            set(findobj(Fmenu,'Tag','Modality'),'Value',ModNum,'UserData',ModNum);
        else
            warning('SPM Menu window not found.');
        end
        
    case 'close'
        Fmenu = spm_figure('FindWin','Menu');
        delete(Fmenu);
        
    otherwise
        error('Unknown action ''%s''.',action);
end


%==========================================================================
%-function Fmenu = create_menu(Vis)
%==========================================================================
function Fmenu = create_menu(Vis)

if ~nargin, Vis = 'on'; end

%-Get size and scalings
%--------------------------------------------------------------------------
PF    = spm_platform('fonts');             %-Font names (for this platform)
FS    = spm('FontSizes');                  %-Scaled font sizes

WS    = spm('WinScale');                   %-Window scaling factors
Rect0 = spm('WinSize','0',1);
RectM = spm('WinSize','Menu',1).*WS;       %-Menu window rectangle: 400x445
Pos   = [Rect0(1) Rect0(2) 0 0] + RectM;

%-Create Menu window
%--------------------------------------------------------------------------
if ~isempty(which('openfig'))
    Fmenu = openfig(fullfile(spm('Dir'),'spm_Menu.fig'),'new','invisible');
else
    fprintf('\nFunction openfig.m is missing. Download it from:\n');
    fprintf(' https://savannah.gnu.org/bugs/download.php?file_id=33499\n');
    Fmenu = [];
end

set(Fmenu,'Name',[spm('Version') ': Menu']);
set(Fmenu,'Units','pixels', 'Position',Pos);
set(Fmenu,'Color',[1 1 1]*.8);

% Fmenu = figure('Tag','Menu',...
% 	'Name',[spm('Version') ': Menu'],...
%     'IntegerHandle','off',...
%     'NumberTitle','off',...
%     'Units','pixels',...
%     'Position',Pos,...
%     'Resize','on',...
%     'Pointer','arrow',...
%     'Color',[1 1 1]*.8,...
%     'MenuBar','none',...
%     'DockControls','off',...
%     'CloseRequestFcn','spm(''Quit'')',...
%     'HandleVisibility','off',...
%     'DefaultUicontrolFontName',PF.helvetica,...
%     'Visible','off');

%-Set SPM colour
%--------------------------------------------------------------------------
set(findobj(Fmenu,'Tag', 'frame'),'BackgroundColor',spm('colour'));
if ismac
    set(findobj(Fmenu,'UserData','LABEL'),'Visible','off','Tag','');
end

%-Set Utils
%--------------------------------------------------------------------------
set(findobj(Fmenu,'Tag', 'Utils'), 'String',{'Utils...',...
    'CD',...
    'PWD',...
    'Run M-file',...
    'Load MAT-file',...
    'Save MAT-file',...
    'Delete files',...
    'Show SPM',...
    'Show MATLAB'});
set(findobj(Fmenu,'Tag', 'Utils'), 'UserData',{...
    ['spm(''FnBanner'',''CD'');' ...
     'cd(spm_select(1,''dir'',''Select new working directory''));' ...
     'spm(''alert"'',{''New working directory:'',[''    '',pwd]},''CD'',1);'],...
    ['spm(''FnBanner'',''PWD'');' ...
     'spm(''alert"'',{''Present working directory:'',[''    '',pwd]},''PWD'',1);'],...
    ['spm(''FnBanner'',''Run M-file'');' ...
     'spm(''Run'');' ...
     'fprintf(''%-40s: %30s\n'',''Completed'',spm(''time''));'],...
    ['spm(''FnBanner'',''Load MAT-file'');' ...
     'load(spm_select(1,''mat'',''Select MAT-file''));'],...
    ['spm(''FnBanner'',''Save MAT-file'');' ...
     'save(spm_input(''Output filename'',1,''s''), spm_get_defaults(''mat.format''));'],...
    ['spm(''FnBanner'',''Delete files'');' ...
     'spm(''Delete'');'],...
    ['spm(''FnBanner'',''Show SPM'');' ...
     'spm(''Show'');'],...
    ['spm(''FnBanner'',''Show Command Window'');' ...
     'commandwindow;']});

%-Set Toolboxes
%--------------------------------------------------------------------------
xTB       = spm('tbs');
if ~isempty(xTB)
    set(findobj(Fmenu,'Tag', 'Toolbox'),'String',{'Toolbox:' xTB.name });
    set(findobj(Fmenu,'Tag', 'Toolbox'),'UserData',xTB);
else
    set(findobj(Fmenu,'Tag', 'Toolbox'),'Visible','off')
end
set(Fmenu,'Visible',Vis);
