function ndx = good_channel(Channel,ChannelFlag,Str);
%GOOD_CHANNEL - Extract channels of a given type
% function ndx = good_channel(Channel,ChannelFlag,Str);
% ndx = good_channel(Channel,ChannelFlag,Str) 
% or 
% ndx = good_channel(Channel,[],Str);
% Return the index to the good channels for desired sensor type
% INPUTS:
%       - CHANNEL is a BrainStorm Channel.mat structure, of importance here is Channel.Type
%       - CHANNELFLAG is vector of -1 (bad), 0(indifferent), and 1(good) data flags of length 
%                     the number of channels in the Channel structure.
%                     If left empty, assumes all channels are GOOD: % Str may contain the name of a channel

%<autobegin> ---------------------- 27-Jun-2005 10:44:23 -----------------------
% ------ Automatically Generated Comments Block Using AUTO_COMMENTS_PRE7 -------
%
% CATEGORY: Forward Modeling
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\good_channel.m  NOTE: Routine calls itself explicitly
%
% At Check-in: $Author: Mosher $  $Revision: 20 $  $Date: 6/27/05 8:59a $
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
%<autoend> ------------------------ 27-Jun-2005 10:44:23 -----------------------


%       - STR is either one of the valid Channel.Type strings:
%             ('MEG','EEG','STIM','OTHER','FUSION','EEG REF', 'MEG REF'), or 'ALL'.
%             STR may also contain the name of a channel 
%             (ex. STR ='MRT32')
% OUPUTS:
%       - ndx, an indexer into Channel(ndx) and Data.F(ndx,:), etc, of good channels
%              Good channels are those that are not bad, i.e. both flags of 0 and 1
% See also GET_CHANNEL

% /---Script Author--------------------------------------\
% |                                                      |
% |  *** John C. Mosher, Ph.D.                           |
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
%---------------------------------------------------------------------------------------------------------------------------
% 2/16/00 JCM changed logic such that good channels are not bad, i.e. both 0 and 1
%  are allowed
% 11/29/00 SB - Updated the help header and allowed empty ChannelFlag vectors
% 11-2001 SB, added the FUSION type to select both the MEG and EEG channels
% SB  18-Nov-2002 Added the 'EEG REF' Channel type to detect the channel number of the EEG reference for a given Channel set
% SB  20-Nov-2002 Added the 'MEG REF' Channel type for fast access to these MEG reference channels during os_meg call and
%                 remove .Gcoef and .i*sens  fields from Channel structure 
% SB 16-Mar-2003  Now return channels with ChannelFlag == 0 too.
% SB 14-Apr-2005  Str may contain the name of a channel
%---------------------------------------------------------------------------------------------------------------------------

Flag = zeros(1,length(Channel)); % one flag per channel
Str = deblank(lower(Str)); % case insensitivity

if isempty(ChannelFlag) % Assume all channels are good
   ChannelFlag = ones(length(Channel),1);
end

   
switch Str
case 'all'
   % get all of the good channels
   ndx = find(ChannelFlag(:)' >= -1);
   return
   
case {'meg','eeg','other','stim','eeg ref','meg ref'}
   
   for i = 1:length(Channel),
      if(strcmp(deblank(lower(Channel(i).Type)),Str)),
         Flag(i) = 1; % it is a match to Str
      end
   end
   
   % now form logical on good channel flag, then accept those 
   %  that are both true in the Flag and True as Good
   ndx = find(all([Flag;(ChannelFlag(:)' >= 0)]));
   return
   
case 'fusion'
    ndxmeg = good_channel(Channel,ChannelFlag,'MEG');
    ndxeeg = good_channel(Channel,ChannelFlag,'EEG');
    ndx = [ndxmeg,ndxeeg];
    
    otherwise
   
    
    try % Str may contain the name of a channel
        ndx = find(strcmpi(Str,{Channel.Name}));
    catch
        errordlg(sprintf('Unknown Type or Name: %s',Str),'Good_Channel','modal');
    end
    
   
end

return
