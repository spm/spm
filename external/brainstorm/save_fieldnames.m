function save_fieldnames(A,fname,append);
%SAVE_FIELDNAMES - Save just the fieldnames of a structure to a mat-file
% function save_fieldnames(A,fname,append);
% The individual fields within A are saved to fname, such that
% A = load(fname) reverts back to the original structure
% when append is '-append', ie the optional string of the Matlab SAVE command,
% fields of A are appended to existing fname

%<autobegin> ---------------------- 27-Jun-2005 10:45:38 -----------------------
% ------ Automatically Generated Comments Block Using AUTO_COMMENTS_PRE7 -------
%
% CATEGORY: Utility - General
%
% At Check-in: $Author: Mosher $  $Revision: 17 $  $Date: 6/27/05 9:00a $
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
%<autoend> ------------------------ 27-Jun-2005 10:45:38 -----------------------

% ----------------------------- Script History ---------------------------------
% John C. Mosher, Ph.D.
% SB 01-Aug-2002 : Added optional 'append' argument
% 19-May-2004 JCM Comments Cleaning
% ----------------------------- Script History ---------------------------------


NAMES = fieldnames(A); % the fieldnames
% extract each fieldname into the workspace
for i = 1:length(NAMES),
   eval(sprintf('%s = getfield(A,''%s'');',NAMES{i},NAMES{i}));
end

if nargin == 2
    save(fname,NAMES{:}) % save the fields to the file
elseif nargin == 3  & strcmp(lower(append),'-append') % Append variable(s)
    save(fname,NAMES{:},'append') % save the fields to the file
end

return
