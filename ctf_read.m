function [ctf] = ctf_read(folder,CHAN,TIME,TRIALS,COEFS);

% ctf_read - Read data from a CTF .ds folder
%
% [ctf] = ctf_read( [folder], [CHAN], [TIME], [TRIALS], [COEFS] );
% 
% eg,
%     ctf = ctf_read; % inputs are all optional
%     ctf = ctf_read('/data/directory/datasetname.ds');
%     ctf = ctf_read(folder,'meg','all','all');
% 
% ctf struct has fields:
%
% ctf.data - all of the data, 'data{1}' will display the first data set
% ctf.header
% ctf.setup
% ctf.sensorNames
% ctf.sensorLocations
% ctf.sensorOrientations
%
% This function calls, in this order:
%
% ctf_read_res4 - to read in header, gain/offset, and sensor information
% ctf_read_meg4 - to read in the data
% 
% INPUTS---------------------------------------------------------------------
% folder:     The .ds directory of the dataset.  If not given, a graphical
%             prompt for the folder appears.
%
% CHAN:       eg. [30:35] - an interval of the desired channels to be read.
%             If CHAN = 'eeg', only eeg channels/sensors
%             If CHAN = 'meg', only meg channels/sensors
%             If CHAN = 'ref', only reference channels/sensors
%             If CHAN = 'other' only the other channels/sensors
%             If CHAN = 'megeeg', only meg and eeg channels/sensors
% 
% TIME:       eg. [0 5] - seconds 0 to 5, the time interval to read.
%             If TIME = 'all', the entire duration of the trial(s) will
%             be read (i.e. TIME = [1:ctf.setup.duration]).
%             
% TRIALS:     eg. TRIALS = n, the nth trial is read.
%             eg. TRIALS = [3,5,8], trials 3, 5, and 8 are read such that,
%                 ctf.data{1} = data for trial 3, 
%                 ctf.data{2} = data for trial 5, and
%                 ctf.data{3} = data for trial 8.
%             eg. TRIALS = [3:7], trials 3 through 7 are read.            
%             eg. TRIALS = 'alltrials', the data for all of the trials are
%                 read (i.e. TRIALS = [1:ctf.setup.duration]).
%             
% COEFS:   an option to read the sensor coefficients, which give the
%          weights for calculation of synthetic 2nd or 3rd order
%          gradiometers.
%          If coefs = 1, read the sensor coefficients
%          If coefs = 0, do not read the sensor coefficients (default)
% 
%
%      <>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> %
%      <                                                      > %  
%      <                    DISCLAIMER:                       > %
%      <                                                      > %
%      < THIS PROGRAM IS INTENDED FOR RESEARCH PURPOSES ONLY. > %
%      < THIS PROGRAM IS IN NO WAY INTENDED FOR CLINICAL OR   > %
%      <                    OFFICIAL USE.                     > %
%      <                                                      > %
%      <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<> %
%



% $Revision: 253 $ $Date: 2004/08/19 03:17:10 $

% Copyright (C) 2004  Darren L. Weber
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Modified: 11/2003, Darren.Weber_at_radiology.ucsf.edu
%                    - modified from NIH code
%                      simply to allocate data into one large struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('folder','var'),
  ctf = ctf_folder;
else
  ctf = ctf_folder(folder);
end

if ~exist('CHAN','var'),   CHAN   = 'all';  end
if ~exist('TIME','var'),   TIME   = 'all';  end
if ~exist('TRIALS','var'), TRIALS = 'all'; end
if ~exist('COEFS','var'),  COEFS  = 0; end

if isempty(CHAN),   CHAN   = 'all';  end
if isempty(TIME),   TIME   = 'all';  end
if isempty(TRIALS), TRIALS = 'all'; end
if isempty(COEFS),  COEFS  = 0; end

ctf = ctf_read_res4(ctf.folder,1,COEFS);

ctf = ctf_read_meg4(ctf.folder,ctf,CHAN,TIME,TRIALS,COEFS);

return
