function nrm = rownorm(A);
%ROWNORM - Calculate the L2 norm of each row of a matrix
% function nrm = rownorm(A);
% calculate the Euclidean norm of each ROW in A, return as
% a column vector with same number of rows as A
%
% See also COLNORM

%<autobegin> ---------------------- 27-Jun-2005 10:45:37 -----------------------
% ------ Automatically Generated Comments Block Using AUTO_COMMENTS_PRE7 -------
%
% CATEGORY: Utility - Numeric
%
% At Check-in: $Author: Mosher $  $Revision: 16 $  $Date: 6/27/05 9:00a $
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
%<autoend> ------------------------ 27-Jun-2005 10:45:37 -----------------------

% ----------------------------- Script History ---------------------------------
% 1994 by John C. Mosher
% May 6, 1994 JCM author
% 19-May-2004 JCM Comments Cleaning
% ----------------------------- Script History ---------------------------------

[m,n] = size(A);

if(n>1),			% multiple columns
  nrm = sqrt(sum([A.*conj(A)],2));
else				% A is col vector
  nrm = abs(A);
end

return
