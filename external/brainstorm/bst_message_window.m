function varargout = bst_message_window(varargin)
%BST_MESSAGE_WINDOW - Application M-file for bst_message_window.fig, with NON CALLBACKS
% function varargout = bst_message_window(varargin)
%    FIG = BST_MESSAGE_WINDOW launch bst_message_window GUI.
%    BST_MESSAGE_WINDOW('callback_name', ...) invoke the named callback.
% NON CALLBACKS:
% bst_message_window(cellstr) appends cell array of strings cellstr to window 
% bst_message_window('append',str) appends string or cell of strings to window
% bst_message_window(str) will also append string, unless it is a valid function call
% bst_message_window('wrap',str) will wrap string or cell of strings to the window
% bst_message_window('unique',str) will wrap string or cell of strings to a unique window,
% bst_message_window('process',str) send str to the ProcessLauncher GUI (see also bst_ProcessLauncher)
%  returning the handle to the window, i.e. it won't go to the message window, but rather
%  its own message window.
% bst_message_window('overwrite',str) overwrites the last line with the new string
% 
% bst_message_window('close') removes the window, otherise the CloseRequestFcn executes
% If window does not exist, 'append' will open it

%<autobegin> ---------------------- 27-Jun-2005 10:43:34 -----------------------
% ------ Automatically Generated Comments Block Using AUTO_COMMENTS_PRE7 -------
%
% CATEGORY: GUI and Related
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\bst_color_scheme.m
%   toolbox\bst_layout.m
%   toolbox\bst_message_window.m  NOTE: Routine calls itself explicitly
%
% Subfunctions in this file, in order of occurrence in file:
%   varargout = process(str,fig);
%   varargout = append(str,fig);
%   varargout = unique(str);
%   varargout = wrap(str,fig);
%   varargout = overwrite(rep_str,fig);
%   varargout = close();
%   varargout = Messages_Callback(h, eventdata, handles, varargin)
%   varargout = clear_selected_Callback(h, eventdata, handles, varargin)
%   varargout = clear_all_Callback(h, eventdata, handles, varargin)
%   varargout = BrainStormMessages_CloseRequestFcn(h, eventdata, handles, varargin)
%   varargout = BrainStormMessages_ResizeFcn(h, eventdata, handles, varargin)
%
% Application data and their calls in this file:
%   'BrainStormMessageWindow'
%   'MAXLength'
%   'TileType'
%   
%   setappdata(0,'BrainStormMessageWindow',fig);
%   setappdata(fig,'MAXLength',200);
%   setappdata(fig,'TileType','M');
%   setappdata(fig,'TileType','T');
%   
%   MAXLength = getappdata(fig,'MAXLength');
%   delete(getappdata(0,'BrainStormMessageWindow'));
%   fig = getappdata(0,'BrainStormMessageWindow');
%   if(~ishandle(getappdata(0,'BrainStormMessageWindow')))
%
% Figure Files opened by this function:
%   mfilename
%
%   Format of strings below: Type:Style:Tag, "String", CallBack Type and Call
%   <automatic> callback is <Tag>_Callback by Matlab default
%
% Callbacks by figure bst_message_window.fig
%   figure::BrainStormMessages "" uses ResizeFcn for <automatic>
%   uicontrol:listbox:Messages "BrainStorm Message Window" uses Callback for <automatic>
%
% At Check-in: $Author: Mosher $  $Revision: 28 $  $Date: 6/27/05 8:59a $
%
% This software is part of BrainStorm Toolbox Version 27-June-2005  
% 
% Principal Investigators and Developers:
% ** Richard M. Leahy, PhD, Signal & Image Processing Institute,
%    University of Southern California, Los Angeles, CA
% ** John C. Mosher, PhD, Biophysics Group,
%    Los Alamos National Laboratory, Los Alamos, NM
% ** Sylvain Baillet, PhD, Cognitive Neuroscience & Brain Imaging Laboratory,
%    CNRS, Hopital de la Salpetriere, Paris, France
% 
% See BrainStorm website at http://neuroimage.usc.edu for further information.
% 
% Copyright (c) 2005 BrainStorm by the University of Southern California
% This software distributed  under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html .
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%<autoend> ------------------------ 27-Jun-2005 10:43:34 -----------------------


% ----------------- Change History ------------------------
% JCM 27 Feb 2002 changed propertyname to BrainStormMessageWindow, handled "catch" for string input
% JCM 3 May 2002 added 'wrap' and 'overwrite' options
% JCM 10 May 2002 changed CloseRequestFcn to allow non-modal closing.
% JCM 15 May 2002 updated 'overwrite' to allow multiple lines of overwrite
% SB  15 Oct 2002 when BrainStorm GUIs are closed (e.g. when using command-line BrainStorm)
%                 bst_message_window uses the basic DISP command to display info in 
%                 Matlab's command window
% JCM 29 Oct 2002 Preferences should not be erased. Test for presence of taskbar, 
%                 using getappdata(0,'BrainStormTaskbar'); Fixed MAXLength use, should be
%                 an appdata in the message figure, also rolls the last MAXLength lines,
%                 rather than clearing. All cases now handle non-gui mode.
%                 Redid message window to be along the bottom of the screen.
% JCM 13 May 2003 Moved message window buttons to new layout_manager, made window
%                 compatible with layout manager
% JCM 09 Jun 2003 Added "unique" function to allow use for a convenient message and information window
% JCM 26 Jan 2005 Felix fixed rep_str typo in overwrite function
% SB  08 Feb 2005 Added the 'process' option: 
%                 -> bst_message_window('process',str): send str to the ProcessLauncher GUI (see also bst_ProcessLauncher)
% JCM 08-Jun-2005 Set fonts in the creation and resize functions to be
%                 fixed 8 point size,after running the bst_color_scheme function.
% ----------------------------------------------------------

% Last Modified by GUIDE v2.5 08-Jun-2005 15:20:15

if nargin == 0  % LAUNCH GUI
   
   
   fig = openfig(mfilename,'reuse');
   
   % Use system color scheme for figure:
   set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));
   
   % Generate a structure of handles to pass to callbacks, and store it. 
   handles = guihandles(fig);
   guidata(fig, handles);
   
   if nargout > 0
      varargout{1} = fig;
   end
   
   % BrainStorm specific code here
   bst_color_scheme(fig);
   set(handles.Messages,'fontname','FixedWidth'); % good for tables
   % lock in font units and size to keep from scaling
   set(handles.Messages,'fontunits','points','fontsize',8);
   
   setappdata(fig,'TileType','M'); % message window
   bst_layout('align',fig);
   
   setappdata(0,'BrainStormMessageWindow',fig); % set the handle in the taskbar application data
   
   % ------ SET MAXIMUM NUMBER OF LINES TO BE DISPLAYED ------------------
   setappdata(fig,'MAXLength',200); % Maximum number of lines to be displayed in the message window.
   % too many lines significantly slows down all BsT processes.
   % JCM changed 29 Oct to roll the last 200 lines, not clear.
   
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
   
   try      
      if (nargout)
         [varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
      else
         feval(varargin{:}); % FEVAL switchyard
      end
   catch
      % string was not a valid function call, treat as a string
      try
         append(varargin); % try appending the string
      catch
         disp(lasterr);
      end
   end
   
elseif iscell(varargin{1}) % input is a cell array
   
   append(varargin{1});
   
end


% ---------------- subfunctions -----------------

% --------------------------------------------------------------------
function varargout = process(str,fig);
% Send messages to ProcessLauncher GUI
processGUI = findobj(0,'type','figure','tag','bst_ProcessLauncher');
HprocessGUI = guihandles(processGUI);
set(HprocessGUI.dProgressReport,'String',str)

% --------------------------------------------------------------------
function varargout = append(str,fig);
% Call as bst_message_window('append',str)
%  alternatively, if string is unrecognized command, switchyard will call this subfunction

if ischar(str),
   str = {str}; % convert to cell array
end

% Test whether user is calling BrainStorm routines from command line function calls or GUIs.

if ~isappdata(0,'BrainStormTaskbar')
   % The taskbar does not exists, assume we are in command line mode
   for i = 1:length(str),
      disp(str{i});
   end
   return % break away
end
   
% else we are in GUI mode

% get the handle to the message window
if(~isappdata(0,'BrainStormMessageWindow')), % missing
   bst_message_window; % create it
end
if(~ishandle(getappdata(0,'BrainStormMessageWindow')))
    bst_message_window; % invalid handle, create it again
end

% append a string to the message window. str is either string or cell of strings.
if(~exist('fig','var')), % caller did not give fig handle (default)
    fig = getappdata(0,'BrainStormMessageWindow');
end

handles = guidata(fig); % handles in that figure

MAXLength = getappdata(fig,'MAXLength');

lenstr = length(str); % number of lines being added

oldstr = get(handles.Messages,'string');
if ischar(oldstr)
   oldstr = {oldstr}; % convert to cell array
end

oldstr(end+[1:lenstr]) = str; % add new strings to end of old string

ndx = [-MAXLength:0]+length(oldstr); % last MAXLength lines

% ndx may be negative, if length shorter than MAXLength
ndx = ndx(ndx > 0); % keep only the positive ones
 
%trim
oldstr = oldstr(ndx);

% select all of the new message for visibility to the user
set(handles.Messages,'string',oldstr,...
   'val',length(oldstr)+[(1-lenstr):0],...
   'listboxtop',length(oldstr)); % put in

figure(fig); % bring to the front

drawnow

if nargout > 0
    varargout{1} = fig;
end



% -------------------------------------------------------------------
function varargout = unique(str);
% wrap text to a unique message window, not the main window
% JCM 9-Jun-2003 creating a variant of the message window

fig = openfig(mfilename,'new'); % new unique filename

% Use system color scheme for figure:
set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

% Generate a structure of handles to pass to callbacks, and store it. 
handles = guihandles(fig);
guidata(fig, handles);

if nargout > 0
    varargout{1} = fig;
end

% BrainStorm specific code here
bst_color_scheme(fig);
set(handles.Messages,'fontname','FixedWidth'); % good for tables
setappdata(fig,'TileType','T'); % tile
bst_layout('align',fig,1,2,1); % default action, can be set in calling routine as well

% customize for a unique message window
set(fig,'Name','Information','Tag','Information','CloseRequestFcn','closereq','ResizeFcn',[]);
overwrite(sprintf('Information Message %s',datestr(now)),fig); % replace the "BrainStorm Message Window" line

% ------ SET MAXIMUM NUMBER OF LINES TO BE DISPLAYED ------------------
setappdata(fig,'MAXLength',200); % Maximum number of lines to be displayed in the message window.
% too many lines significantly slows down all BsT processes.
% JCM changed 29 Oct to roll the last 200 lines, not clear.

wrap(str,fig); % now wrap the text to just this unique fig



% --------------------------------------------------------------------
function varargout = wrap(str,fig);
% called as mfilename('wrap',str), wrap string to the window

if ischar(str),
   str = {str}; % convert to cell for consistent handling
end

if ~isappdata(0,'BrainStormTaskbar')
   % The taskbar does not exists, assume we are in command line mode
   for i = 1:length(str),
      disp(str{i});
   end
   return % break away
end

% else we are in GUI Mode
% get the handle to the message window

if(~isappdata(0,'BrainStormMessageWindow')), % missing
   bst_message_window; % create it
end
if(~exist('fig','var')), % user did not give as input to this function (default)
    fig = getappdata(0,'BrainStormMessageWindow');
end

handles = guidata(fig); % handles in that figure

if(0), % deprecated, user may send cell array of strings
   % first convert str into one long message string, with proper spaces
   msg_str = [];
   for i = 1:length(str),
      msg_str = [msg_str str{i}];
      if(~strcmp(msg_str(end),' ')),
         msg_str(end+1) = ' ';
      end
   end
   msg_str(end) = []; % remove last space
else
   msg_str = str; % don't alter
end

% wrap the text to the message window
% textwrap in R12.1 is a little too wide
mesPos = get(handles.Messages,'position');
mesWidth = mesPos(3); % current width
set(handles.Messages,'position',[mesPos(1:2) mesWidth*.95 mesPos(4)]); % slightly narrower
outstring = textwrap(handles.Messages,msg_str); % wrap to columns
set(handles.Messages,'position',mesPos); % original size

append(outstring,fig);

if nargout > 0
    varargout{1} = fig;
end



% --------------------------------------------------------------------
function varargout = overwrite(rep_str,fig);
% Overwrite the last lines of the message window with cell array of rep_str
%  if rep_str is a string, overwrite just the last line.
% Useful for updating a "Processing . . ." with a "Done"
%  or a pseudo waitbar,
% bst_message_window('Wait: 1 of 10')
% for i = 2:10,bst_message_window('overwrite',sprintf('Wait: %.0f of 10',i)),<process>,end

if(ischar(rep_str)),
   rep_str = {rep_str}; % make cell for consistent handling
end


if ~isappdata(0,'BrainStormTaskbar')
   % The taskbar does not exists, assume we are in command line mode
   for i = 1:length(rep_str),
      disp(rep_str{i});
   end
   return % break away
end


NumLines = length(rep_str); % how many lines to replace

% get the handle to the message window
if(~isappdata(0,'BrainStormMessageWindow')), % missing
   bst_message_window; % create it
end

if(~exist('fig','var')), % user did not give as input to this function (default)
    fig = getappdata(0,'BrainStormMessageWindow');
end
handles = guidata(fig); % handles in that figure

str = get(handles.Messages,'string'); %existing strings

if(ischar(str)),
   str = {str}; % make cell for consistent handling
end

if(length(str) >= NumLines),
   str([(1-NumLines):0]+end) = rep_str;
else
   % there are not enough lines in the display, just append
   str(end+[1:NumLines]) = rep_str; 
end

set(handles.Messages,...
   'string',str,'listboxtop',...
   length(str),'val',length(str)+[(1-NumLines):0]); % put in, with new strings clearly showing

if nargout > 0
    varargout{1} = fig;
end




% --------------------------------------------------------------------
function varargout = close(); 
% close the message window. Otherwise, the CloseRequestFcn will execute
delete(getappdata(0,'BrainStormMessageWindow')); % close the window
rmappdata(0,'BrainStormMessageWindow'); % clear the application data


% -------------- Callback routines -----------------------

% --------------------------------------------------------------------
function varargout = Messages_Callback(h, eventdata, handles, varargin)
% double clicking will remove a line. Otherwise, no effect
% doulbe right clicking will remove lots of highlighted lines

status = get(handles.BrainStormMessages,'SelectionType');
switch status
case 'normal' % do nothing
case 'open' % user double clicked
   clear_selected_Callback(h,eventdata,handles,varargin);
end





% --------------------------------------------------------------------
function varargout = clear_selected_Callback(h, eventdata, handles, varargin)
% clear the selected messages, keep the view stable

val = get(handles.Messages,'val'); % vector of lines to delete
str = get(handles.Messages,'string'); %existing strings
listboxtop = get(handles.Messages,'listboxtop'); %what is showing at the top of window

str(val) = []; % remove the existing strings
val = val(1); % set to last string deleted
val = min(val,length(str)); % in acceptable range

listboxtop = min(listboxtop,val); % in acceptable range
listboxtop = min(listboxtop,length(str)); % in acceptable range

if(length(str) == 0),
   set(handles.Messages,'val',1,'ListBoxTop',1,'string',{'BrainStorm Message Window'})
else
   set(handles.Messages,'val',val,'string',str,'listboxtop',listboxtop);
end





% --------------------------------------------------------------------
function varargout = clear_all_Callback(h, eventdata, handles, varargin)

% no longer a button for this 13 May 2003
% can instead highlight lots of lines, then double right click to remove them
ButtonName = questdlg('Clear all lines?','Message Window','Yes','No','No');

switch ButtonName
case 'Yes'
   set(handles.Messages,'val',1,'ListBoxTop',1,'string',{'BrainStorm Message Window'})
end





% --------------------------------------------------------------------
function varargout = BrainStormMessages_CloseRequestFcn(h, eventdata, handles, varargin)

% since this window now automatically opens when needed, this forced closure
%  is no longer really necessary. JCM 13-May-2002

if(0), % deprecated code
   str = {'','This BrainStorm message window should stay open.',''}; 
   bst_message_window('append',str);
   
   ButtonName = questdlg('Close the message window?','Message Window','Yes','No','No');
   
   switch ButtonName
   case 'Yes'
      bst_message_window('close');
   end
else % newer code, non-modal closing
   bst_message_window('close');
   disp(' ')
   disp('Generally, the BrainStorm Message Window should remain open.')
   disp('Message window closed.')
end



% --------------------------------------------------------------------
function varargout = BrainStormMessages_ResizeFcn(h, eventdata, handles, varargin)

% want to keep the fonts at the same size

bst_color_scheme(handles.BrainStormMessages);

set(handles.Messages,'fontunits','points');
set(handles.Messages,'fontsize',8);
