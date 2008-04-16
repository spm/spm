function d = dist_sph(vec,sensloc)
%DIST_SPH - FMINS distance function used to minimize the fit to a sphere
% function d = dist_sph(vec,sensloc)
% Given center and list of sensor points
%  find the average distance from the center to these points

%<autobegin> ---------------------- 27-Jun-2005 10:44:11 -----------------------
% ------ Automatically Generated Comments Block Using AUTO_COMMENTS_PRE7 -------
%
% CATEGORY: Utility - Numeric
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\norlig.m
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
%<autoend> ------------------------ 27-Jun-2005 10:44:11 -----------------------

% History -----------------------------------------
% JCM mid-1990s  creation
% SB late 90s  update
% JCM 08-Sep-2003  comments
% -------------------------------------------------

R = vec(end);
center = vec(1:end-1);
nmes = size(sensloc,1);
d = mean(abs(norlig(sensloc - ones(nmes,1)*center) - R)); % Average distance between the center if mass and the electrodes
