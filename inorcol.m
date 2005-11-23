function [normeM] = inorcol(Mat)
%INORCOL - Compute the (pseudo)inverse of the column norms of matrix Mat
% function [normeM] = inorcol(Mat)
% normeM is a sparse diagonal matrix whose diagonal elements 
%  are the inverse of the corresponding column norm of matrix Mat.
%  If a column is zero, then its inverse is set to zero as well.

%<autobegin> ---------------------- 27-Jun-2005 10:44:53 -----------------------
% ------ Automatically Generated Comments Block Using AUTO_COMMENTS_PRE7 -------
%
% CATEGORY: Utility - Numeric
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
%<autoend> ------------------------ 27-Jun-2005 10:44:53 -----------------------

% ----------------------- History -----------------------------------
% JCM 19-May-2003  Handling error case of zero column norm
% ----------------------- History -----------------------------------


cn = sqrt(sum(Mat.*Mat,1)'); % sum only the rows into a column vector
Zero = max(cn)*eps; % relative concept of zero
ndx = cn > Zero;
cn(ndx) = 1 ./ cn(ndx); % zero values are kept zero for a pseudoinverse

normeM = spdiags(cn,0,length(cn),length(cn));

if(0) % old code
    tmp = speye(size(Mat,2));
    normeM = 1./sqrt(sum(Mat.*Mat))';
    normeM = spdiags(normeM,0,tmp);
end
