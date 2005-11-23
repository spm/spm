function Users = get_user_directory;
%GET_USER_DIRECTORY - Get the user's root directory information
% function Users = get_user_directory;
% User = get_user_directory;
% Accesses the Preferences to get the current Users structure and
%  CurrentData, as a subfield fo Users, e.g. Users.CurrentData.
%
% Every attempt to get the current user database should call this mfile, which will in 
%  in turn call database_manager only if there is a problem.
%
% LOGIC: if the user has an active BrainStorm session running, then the structure
%  will be found in the BrainStormTaskbar as application data BrainStormDataBase.
% If there is something wrong with that information, send the user to the data manager
%
% If the BrainStormTaskbar is not running, then we may be in batch mode, so get
%  the information from the BrainStorm preferences. Send error messages to the
%  standard output.
%
% Returns the structure with fields:
%  Comment, string commenting on the location of the files
%  STUDIES, string giving the root location of all studies
%  SUBJECTS, string giving the root for subjects.
%  FILELIST, structure, result of browse_study_folder(SUBJECTS)
%  CurrentData with fields and examples:
%      StudyFile: 'jcm_spont08-31.ds\jcm_brainstormstudy.mat'
%    SubjectFile: 'jcm_mri\jcm_brainstormsubject.mat'
%
% Idea is that we then call load(fullfile(User.STUDIES,<filename>));
%  i.e. all filenames are referential to a root.

%<autobegin> ---------------------- 27-Jun-2005 10:44:22 -----------------------
% ------ Automatically Generated Comments Block Using AUTO_COMMENTS_PRE7 -------
%
% CATEGORY: Utility - General
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\browse_study_folder.m
%   toolbox\bst_message_window.m
%   toolbox\data_manager.m
%
% Group : Preference data and their calls in this file:
%   CurrentData = getpref('BrainStorm','CurrentData',[]);
%   CurrentData = getpref('BrainStorm','UserCurrentData');
%   UserDB = getpref('BrainStorm','UserDataBase');
%   iUser = getpref('BrainStorm','iUserDataBase');
%
% Application data and their calls in this file:
%   TASKBAR = getappdata(0,'BrainStormTaskbar');
%   Users = getappdata(TASKBAR,'BrainStormDataBase');
%
% At Check-in: $Author: Mosher $  $Revision: 24 $  $Date: 6/27/05 8:59a $
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
%<autoend> ------------------------ 27-Jun-2005 10:44:22 -----------------------


% /---Script Author--------------------------------------\
% | *** John C. Mosher, Ph.D.                            |
% |  Biophysics Group                                    |
% |                                                      |
% | *** Sylvain Baillet, Ph.D.                           |
% | Cognitive Neuroscience & Brain Imaging Laboratory    |
% | CNRS UPR640 - LENA                                   | 
% | Hopital de la Salpetriere, Paris, France             |
% | sylvain.baillet@chups.jussieu.fr                     |
% \------------------------------------------------------/
%  
% Date of creation: January 1999
% Date of modification: September 2001                  
%----------------------------------------------------------------------------------
% JCM 11/9/99     appended filesep if missing
% JCM 10/25/00    changed to use get_user_directory
% SB 11/13/00     Added an error message to force the user to move Users.mat into 
%                 custom if necessary.
% SB 09/03/01     USer guidata to look for User information
% JCM 22-Mar-2002 allowing for get_user_directory to load from bst_prev_session
% JCM 28-May-2002 database_manager has replaced initialize_users, database_manager now
%                 uses getprefs rather than users.mat. Altered help block above to
%                 reflect new changes. Used try, catch to handle errors
% SB  03-Sep-2002 Return [] when database is not loaded properly 
% JCM 30-Oct-2002 Strict enforcement of Preferences only, all other use disabled.
%                 Since startup sets a default if needed, errors should be minimal.
% JCM 20-Aug-2003 empty CurrentData default if necessary
% JCM  5-Sep-2003 Added FILELIST generation if it is missing from the database
% ----------------------------------------------------------------------------------

try 
   UserDB = getpref('BrainStorm','UserDataBase');
   iUser = getpref('BrainStorm','iUserDataBase'); % number of the last one used
   Users = UserDB(iUser);
   
   if(~isfield(Users,'FILELIST')),
      % FILELIST is missing from the user database structure
      Users = setfield(Users,'FILELIST',browse_study_folder(Users.STUDIES));
   end
   
   % the current data is stored in preferences separate from the database
   CurrentData = getpref('BrainStorm','CurrentData',[]); % default blank if not there
   Users = setfield(Users,'CurrentData',CurrentData); % May BST MMII format
   
catch
   Users = [];
end

return


% everything else below disabled by the above return statement

try
   
   if(isappdata(0,'BrainStormTaskbar')), % a BrainStorm Taskbar is running
      try
         
         TASKBAR = getappdata(0,'BrainStormTaskbar');
         Users = getappdata(TASKBAR,'BrainStormDataBase');
         % NOTE: May BST MMII returns Users.CurrentData as one of the fields, which is
         %  not in the regular structure.
         
      catch
         
         bst_message_window('wrap',{'Something failed in retrieving your latest database.',...
               'Please edit your database and/or select valid data'});
         data_manager; % non-modal call
         return
         
      end
      
   elseif(ispref('BrainStorm','UserDataBase')),
      
      % there is no Taskbar, we may be in batch mode
      % get the latest database information from the Mathworks preferences
      
      try
         
         UserDB = getpref('BrainStorm','UserDataBase');
         iUser = getpref('BrainStorm','iUserDataBase'); % number of the last one used
         Users = UserDB(iUser);
         % the current data is stored in preferences separate from the database
         CurrentData = getpref('BrainStorm','UserCurrentData');
         Users = setfield(Users,'CurrentData',CurrentData); % May BST MMII format
         
      catch
         
         % use fprintf, since brainstorm is apparently not running
         %  and we may be writing standard output out to a file
         fprintf(['Error in loading your database preferences.\n',...
               'BrainStorm is not running, nor is there a ',...
               'valid database structure in the preferences. \n',...
               'Please restart BrainStorm and use the ',...
               'Data Manager to build a proper database.\n']);
         Users = [];
         return
         
      end
      
      
   else 
      
      % The taskbar isn't running, and there are no preferences
      bst_message_window('wrap',{'BrainStorm is not running, nor is there a ',...
            'valid database structure. Please restart BrainStorm and use the ',...
            'Data Manager to build a proper database.'});
      
      return
      
   end
   
   
catch
   
   % if all else fails, call the data manager
   bst_message_window('wrap',{'Something failed in retrieving your latest database.',...
         'Please edit your database and/or select valid data'});
   data_manager; % non-modal call
   return
   
end


% We have a single Users structure, clean up

% some of our commands anticipate at separator at the end, give one:
if(~strcmp(Users.STUDIES(end),filesep)),
   Users.STUDIES(end+1)=filesep;
end

if(~strcmp(Users.SUBJECTS(end),filesep)),
   Users.SUBJECTS(end+1)=filesep;
end

return



% Deprecated code below, not in use as of 3-Sep-2001 (apparently)
if(0), % deprecated code section, marked more clearly 20-Oct-2002 JCM
   
   if(0), %old SB code
      hf = findobj(0,'Tag','TASKBAR'); % where is the taskbar
      Users = guidata(hf);
   end
   
   
   
   if isempty(hf) % User does not use the taskbar figure (may be running bacth processing ?)
      global Users
      [brainstorm_home,tmp] = which('startup'); % find the brainstorm home folder in the user path
      if isempty(brainstorm_home) % Path does not include the BrainStorm folder containing the 'startup' file
         Users.SUBJECTS = pwd;
         Users.STUDIES = pwd;           
      else 
         User.SUBJECTS = brainstorm_home;
         User.STUDIES = brainstorm_home;
      end
   end
   
   default_custom = what('Custom'); % find the custom folder in the users path
   fname = fullfile(default_custom.path,'users.mat'); % the fully qualified name
   if(~exist(fname)),
      %SB 11/13/00 : Error message is necessary for a proper handling of the situation by the user
      errordlg(sprintf('It is necessary to move the Users.mat file in %s first - Please restart BrainStorm when Users.mat is moved.',default_custom.path),'New BrainStorm Version');
      delete(findobj(0,'Tag','TASKBAR'))
      Users = [];   
      return
   else
      load(fullfile(default_custom.path,'users.mat'),'-mat') %  loads the custom folder file 'users.mat'
   end
   
   % the local directory now has the structure Users, comprising .Comment, .STUDIES, and .SUBJECTS
   
   % trim to just the FIRST one
   
   Users = Users(1);
   
   % some of our commands anticipate at separator at the end, give one:
   if(~strcmp(Users.STUDIES(end),filesep)),
      Users.STUDIES(end+1)=filesep;
   end
   
   if(~strcmp(Users.SUBJECTS(end),filesep)),
      Users.SUBJECTS(end+1)=filesep;
   end
   
   return
   
end % deprecated code section
