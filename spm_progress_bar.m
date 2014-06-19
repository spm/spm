function spm_progress_bar(action,varargin)
% Display a 'Progress Bar' in the 'Interactive' window
% FORMAT spm_progress_bar('Init',height,xlabel,ylabel,flgs)
% Initialise the bar in the 'Interactive' window.
% If flgs contains a 't', then use tex interpreter for labels.
%
% FORMAT spm_progress_bar('Set',value)
% Set the height of the bar itself.
%
% FORMAT spm_progress_bar('Set','xlabel',xlabel)
% FORMAT spm_progress_bar('Set','ylabel',ylabel)
% Set the progress bar labels.
%
% FORMAT spm_progress_bar('Set','height',height)
% Set the height of the progress bar.
%
% FORMAT spm_progress_bar('Clear')
% Clear the 'Interactive' window.
%__________________________________________________________________________
% Copyright (C) 1996-2013 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_progress_bar.m 6057 2014-06-19 11:49:38Z guillaume $

if ~nargin, action = 'Init'; end

% Find the Interactive window and exit if not
%--------------------------------------------------------------------------
Finter = spm_figure('FindWin','Interactive');
if isempty(Finter), return; end

switch lower(action)
    % Initialise
    %----------------------------------------------------------------------
    case 'init'
        if nargin > 1, arg1 = varargin{1}; else arg1 = 1;           end
        if nargin > 2, arg2 = varargin{2}; else arg2 = 'Computing'; end
        if nargin > 3, arg3 = varargin{3}; else arg3 = '';          end
        if nargin > 4, arg4 = varargin{4}; else arg4 = ' ';         end
        if any(arg4 == 't'), interp = 'tex'; else interp = 'none';  end
        pb = struct('pointer',get(Finter,'Pointer'),...
                    'name'   ,get(Finter,'Name'),...
                    'buffer', get(Finter,'DoubleBuffer'));
        spm_progress_bar('Clear');
        set(Finter,'Pointer','watch');
        set(Finter,'Name',pb.name);
        set(Finter,'DoubleBuffer','on');
        pb.ax = axes('Position', [0.45 0.2 0.05 0.6],...
                     'XTick',    [],...
                     'Xlim',     [0 1],...
                     'Ylim',     [0 max([arg1 eps])],...
                     'Box',      'on',...
                     'Parent',   Finter);
        try, set(pb.ax,'ClippingStyle','rectangle'); end
        lab = get(pb.ax,'Xlabel');
        set(lab,'string',arg2,'FontSize',10,'Interpreter',interp);
        lab = get(pb.ax,'Ylabel');
        set(lab,'string',arg3,'FontSize',10,'Interpreter',interp);
        lab = get(pb.ax,'Title');
        set(lab,'string','0% Complete','Interpreter',interp);
        t = clock;
        str = sprintf('Began %2.0f:%02.0f:%02.0f',t(4),t(5),t(6));
        text(2,arg1/2,0,str,'FontSize',10,'Parent',pb.ax);
        l = line('Xdata',     [0.5 0.5],...
                 'Ydata',     [0 0],...
                 'LineWidth', 8,...
                 'Color',     [1 0 0],...
                 'Tag',       'ProgressBar',...
                 'Parent',    pb.ax);
        set(l,'UserData',pb);
        drawnow;
        
    % Set
    %----------------------------------------------------------------------
    case 'set'
        if nargin == 1, value = 0; else  value = varargin{1}; end
        br = findobj(Finter,'Tag','ProgressBar');
        if ~isempty(br)
            pb = get(br,'UserData');
            if ischar(value)
                if nargin == 2, str = ''; else str = varargin{2}; end
                switch lower(value)
                    case {'xlabel','ylabel'}
                        set(get(pb.ax,value),'String',str);
                    case 'height'
                        set(pb.ax,'YLim',[0 max([str eps])]);
                    otherwise
                        error('Unknown action.');
                end
            else
                set(br,'Ydata',[0 value]);
                lim = get(get(br,'Parent'),'Ylim');lim=lim(2);
                lab = get(pb.ax,'Title');
                set(lab,'string',sprintf('%.0f%% Complete',100*value/lim));
            end
            drawnow;
        end
    
    % Clear
    %----------------------------------------------------------------------
    case 'clear'
        pb = get(findobj(Finter,'Tag','ProgressBar'),'UserData');
        spm_figure('Clear',Finter);
        if isstruct(pb)
            set(Finter,'Pointer',     pb.pointer);
            set(Finter,'Name',        pb.name);
            set(Finter,'DoubleBuffer',pb.buffer);
        end
        drawnow;
    
    % Error
    %----------------------------------------------------------------------
    otherwise
        error('Unknown action string');
end
