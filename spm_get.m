function [P] = spm_get(n,Suffix,Prompt,Prefix,NewWDir,CmdLine)
% User interface : filename selection
% FORMAT [P] = spm_get(n,Suffix,Prompt,Prefix);
% n          - 0 - returns with empty P
%            - Positive integer - prompts for n filepaths
%            - +Inf - prompts for filepaths until "Done" pressed
%            - Negative integer - prompts for n directories
%            - -Inf - prompts for directories until "Done" pressed
% Suffix     - Filename filter {i.e. *Suffix} 
%            - numeric 0 for previous filtering
% Prompt     - Prompt string
% Prefix     - Prefix Filename filter {i.e. Prefix*}
% NewWDir    - New working directory
% CmdLine    - CmdLine override switch, 0 for GUI, 1 for Command line.
%
% P      - string matrix of full pathnames {one string per row}
%____________________________________________________________________________
%
% spm_get allows interactive selection of filepaths/directory names
% from disk. Pass positive n for filepaths, negative n for directories.
% +Inf/-Inf for selection of filepaths/directories until "Done" is pressed.
%
% if global CMDLINE exists and is true, the command line is used. The
% CmdLine parameter overrides this choice.
%
% Otherwise, selection of files is via GUI. The select files window, enables
% the interactive directory display and selection of filenames.
% Clicking on a filename (black) adds it to the list of those selected,
% the filename is numbered and changes color. The last numbered filename
% can be unselected by clicking it again. If a specific number of
% filenames have been requested (n), then only this number can be
% selected. Click "Done" when you're finished. Clicking on a directory
% name (red) changes to that directory. The file/directory window is
% scrolled up and down with the buttons provided. The directory window is
% editable, and a pull down menu of previous directories is maintained.
% All files in the current window can be selected *in the order in which
% they appear* with the "All" button, and the "Reset" button clears the
% current list. The "Edit" button enables interactive editing of the
% filename list, and the "Keybd" button invokes a command line
% interface.
%
% The Filter string for filenames is built as Prefix*Suffix, and appears in
% a single editable window. Subsequent calls to spm_get_cb with the same 
% Filter string do not change the value of the Filter string, allowing
% the users edits of the Filer string to persist.
%
% The 'Select Files Window', created by the first call to spm_get, is
% 'Tag'ged 'SelFileWin'. This window is hidden at the end of the spm_get
% transaction, and is used by subsequent calls. This saves time, and
% also permits the storage of previously visited directories from
% transaction to transaction.
%
% CallBacks are handled by spm_get_cb, which requires spm_list_files.c
% (SPM94 routines spm_ls and spm_all are now obsolete.)
%
%__________________________________________________________________________
% %W% Andrew Holmes, Karl Friston %E%


% Version history
% - Karl Friston  - V1 - ??/94
% - Andrew Holmes - V2 - 02/95 - Used Tag objects, no globals.
%                                Directory history, editing, command line &
%                                and unselection of filenames. Doesn't cd!
%                                Negative n for directory selection.

%-Condition arguments
%-----------------------------------------------------------------------------
if nargin<6 CmdLine=[]; end
if nargin<5 NewWDir=''; end
if nargin<4 Prefix=''; end
if nargin<3 Prompt='Select files...'; end
if nargin<2 Suffix=0; end
if nargin<1 n=+Inf; end

if (n==0), P=[]; return, end

if isempty(CmdLine)
	global CMDLINE
	if ~isempty(CMDLINE), CmdLine = CMDLINE; else, CmdLine=0; end
end

%-NB: spm_get_cb uses a single filter string, with 0 for previous filtering
if isstr(Suffix)
	if Suffix(1)~='*', Suffix = ['*',Suffix]; end
	Filter = [Prefix,Suffix];
else
	Filter=0;
end

%-Computation
%=============================================================================

if CmdLine

	%-Command line version, if global CMDLINE is true, & not overridden
	%===================================================================

	P=spm_get_cb('CmdLine',n,Prompt);

else

	%-GUI version
	%===================================================================

	%-Save current figure number
	%-------------------------------------------------------------------
	fig = get(0,'Children'); if ~isempty(fig), fig=gcf; end

	%-Set up FileSelWin
	%-------------------------------------------------------------------
	F = spm_get_cb('Initialise','on',n,Prompt,Filter,NewWDir);
	
	%-Wait until filenames have been selected
	%-------------------------------------------------------------------
	hDone = findobj(F,'Tag','Done');
	while ~get(hDone,'UserData'), pause(0.01), end
	
	%-Recover P
	%-------------------------------------------------------------------
	P=get(findobj(F,'Tag','P'),'UserData');
	
	%-Reset and hide SelFileWin
	%-------------------------------------------------------------------
	spm_get_cb('Initialise','off');

	%-Return focus to previous figure (if any)
	%-------------------------------------------------------------------
	% if ~isempty(fig), figure(fig); end

end % (if CmdLine)

%-Log the transaction
%---------------------------------------------------------------------------
if exist('spm_log')==2
	spm_log(['spm_get : ',Prompt,' :'],P); end
