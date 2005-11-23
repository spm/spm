function [nrm,Anrm] = colnorm(A);
%COLNORM - Calculate L2 norm for each column of a matrix and normalize
% function [nrm,Anrm] = colnorm(A);
% calculate the Euclidean norm of each COLUMN in A, return as
% a row vector with same number of cols as A.
% Optionally, return A with each column now normalized

%<autobegin> ---------------------- 27-Jun-2005 10:43:52 -----------------------
% ------ Automatically Generated Comments Block Using AUTO_COMMENTS_PRE7 -------
%
% CATEGORY: Utility - Numeric
%
% At Check-in: $Author: Mosher $  $Revision: 17 $  $Date: 6/27/05 8:59a $
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
%<autoend> ------------------------ 27-Jun-2005 10:43:52 -----------------------

% ----------------------------- Script History ---------------------------------
% 1994 by John C. Mosher
% May 6, 1994 JCM author 
% Nov 19, 1996 added Anrm option
% June 13, 1999 JCM optimized the nrm division into A
% 19-May-2004 JCM Comments Cleaning
% ----------------------------- Script History ---------------------------------

[m,n] = size(A);

if(m>1),			% multiple rows
  nrm = sqrt(sum([A.*conj(A)]));
else				% A is row vector
  nrm = abs(A);			% just return mag of each column
end

if(nargout > 1),
  ndx = find(nrm>0);		% any zero norm?
  Anrm = zeros(size(A));
  % normalize any non-zero columns
  Anrm(:,ndx) = A(:,ndx) ./ nrm(ones(m,1),ndx);
end

return
