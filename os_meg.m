function G = os_meg(L,Channel,Param,Order,imegsens,irefsens,Verbose);
%OS_MEG - Calculate the (overlapping) sphere models for MEG
% function G = os_meg(L,Channel,Param,Order,imegsens,irefsens,Verbose);
% function G = os_meg(L,Channel,Param,Order,imegsens,irefsens);
% function G = os_meg(L,Channel,Param,Order);
% Modified for CME, not MME, as Order = 1.
% Calculate the magnetic field, spherical head, arbitrary orientation
%
% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% L : a 3 x nL array, each column a source location (x y z coordinates); nL sources
% Channel is a BrainStorm channel structure
% Param an array of structures:
% For every channel i:
% Param(i).Center = a vector of the x, y, z locations for the sphere model 
% (assume the same center for every sphere for the classical spherical head model);
%
% Param(i).Radii = a vector containing the radius in meters of the concentric spheres  ;
% Can be a scalar for the single-sphere head model
%
% Param(i).EEGType = [] % Leave it empty for MEG;       
%
% Order: Defines the source order for which to compute the forward problem:
%  -1 current dipole
%   0 focal(magnetic) dipole
%   1 1st order current multipole
%
% imegsens is the index to the MEG sensors in the Channel information
% irefsens is the index to the MEG reference sensors in the Channel
% if imegsens (irefsens) is not given, then routine (expensively)
% searches the Channel structure for 'MEG' ('MEG REF') values
%
% Verbose : toggle verbose mode (1 is default)
%
% OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% G is the gain matrix: each column is the forward field of each source

%<autobegin> ---------------------- 27-Jun-2005 10:45:21 -----------------------
% ------ Automatically Generated Comments Block Using AUTO_COMMENTS_PRE7 -------
%
% CATEGORY: Forward Modeling
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\blk_diag.m
%   toolbox\blk_lex.m
%   toolbox\bst_message_window.m
%   toolbox\colnorm.m
%   toolbox\cross_mat.m
%   toolbox\good_channel.m
%   toolbox\inorcol.m
%   toolbox\makeuswait.m
%   toolbox\mby3check.m
%   toolbox\norcol.m
%   toolbox\vec.m
%
% Subfunctions in this file, in order of occurrence in file:
%   c = cross(a,b);
%   k = kronmat(a,b);
%   G = sarvas(L,P,Order);
%   G = sarvas_dipole(L,P,Order);
%   D = sarvas_partial(L,P);
%
% At Check-in: $Author: Mosher $  $Revision: 35 $  $Date: 6/27/05 9:00a $
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
%<autoend> ------------------------ 27-Jun-2005 10:45:21 -----------------------

% NOT OPTIMIZED FOR BRAINSTORM STRUCTURE.
%  But it's a whole lot better now (October 99 version).

% John C. Mosher, Ph.D.

% ------------------------ HISTORY -------------------------------------------
% 29-Oct-1999 JCM tried to speed up cross products, which are a substantial portion of 
%                 the total calculation time, for a lot of dipoles.
%                 Also added the waitbar feature for big calculations.
%
% 28-Nov-2000 SB  Modified for updated 3rd-order gradient correction of CTF data  
% 24-Jan-2002 JCM in Paris, modified code to run old partial gradient calculation
%                 in this code
% 20-Feb-2002 JCM fixed missing translocation of the partial gradient location
%                 information. Put center shifting into the partial_gradient 
%                 code itself. Major change to gradient interface to handle
%                 head center in different location.
% 20-Feb-2002 SB  Fixed different version to handle case of m x 3 L matrix,
%                 force L to be 3 x N or return error.
% 12-Mar-2002 JCM Force L to be 3 x N or return error.
% 19-Jun-2002 JCM Realized version differences between toolbox and developers
%                 toolbox, carefully merged two version together into toolbox,
%                 updated comments, history, switched waitbar to message window
% 02-Jul-2002  SB Fixed display using the 'overwrite' option of bst_message_window
% ...........     Added progression messages during computation of large CME gain matrices
% 03-Sep-2002  SB Fixed bug when Channel contains magnetometers and when head center is at origin ([0 0 0])
% 20-Nov-2002  SB Updated computation with CTF 3rd order gradient correction
%                 Now must pass the entire Channel structure to OS_MEG
%                 MEG and MEG REF channels are extracted using GOOD_CHANNEL
%                 Lightly altered display  
% 21-Oct-2003 JCM added optional imegsens and irefsens indexed inputs
% 09-Mar-2004 SB  added Verbose argument
% -----------------------------------------------------------------------------

% Which verbose mode ?
if 0  % deprecated code - get a warning in Matlab 6.5.0 because 'Verbose' is an argument and has to be declared as GLOBAL before first use
    global Verbose % Pass Verbose to subfunctions
end


if nargin == 5
    Verbose = imegsens; 
    clear imegsens
elseif nargin < 7 ^ ~exist('Verbose','var')
    Verbose = 1;% Default
end


% Indices of MEG and MEG REF channels in Channel structure array:
if ~exist('imegsens','var') | isempty(imegsens),
   imegsens = good_channel(Channel,[],'MEG');
end

if isempty(imegsens)
    error('No MEG channels available')
end

if ~exist('irefsens','var'),
   irefsens = good_channel(Channel,[],'MEG REF');
end

if ~isempty(irefsens)
    ChannelRef = Channel(irefsens); % Fill out a specific channel structure with MEG REF channels only
    [ChannelRef(:).Type] = deal('MEG'); % Ref channels will be treated as regular MEG sensors when Sarvas is called to compute their forward fields
    refFlag = 1;
else
    refFlag = 0;
end
Channel = Channel(imegsens); % Keep Channel as a MEG-channel only channel set. 
Param = Param(imegsens);

% ----------------------------

NumCoils = size(Channel(1).Loc,2); % number of coils, assumed same for all channels
NumSensors = length(Channel); % how many channels

if(size([Channel.Loc],2) ~= NumSensors*NumCoils),
   errordlg({'Sorry, OS_MEG not equipped to handle different number of coils'});
   G = [];
   return
end

% load up the old parameter array
% P.sensor is 3 x nR,each column a sensor location
% P.orient is 3 x nR, the sensor orientation
% P.center is 3 x nR, the sphere center for each sensor

[P(1:NumCoils)] = deal(struct('sensor',zeros(3,NumSensors),...
   'orient',zeros(3,NumSensors),...
   'center',zeros(3,NumSensors),'weight',[]));
AllLocs = [Channel.Loc]; % remap all Locations
AllLocs = reshape(AllLocs,NumCoils*3,size(AllLocs,2)/NumCoils);
AllOrient = [Channel.Orient];
AllOrient = AllOrient*inorcol(AllOrient);
AllOrient= reshape(AllOrient,NumCoils*3,size(AllOrient,2)/NumCoils);
AllWeight = [Channel.Weight];
AllWeight = reshape(AllWeight(:),NumCoils,length(AllWeight(:))/NumCoils);
% -- modified from original version: JM 30/10/05 --
AllCenter = [Param.Center];
%AllCenter = reshape(AllCenter,NumCoils*3,size(AllLocs,2)/NumCoils);
% Bug fix by Rik Henson 6/6/07
AllCenter = reshape(AllCenter,NumCoils*3,size(AllCenter,2)/NumCoils);

for j = 1:NumCoils,
   P(j).sensor = AllLocs([-2:0]+j*3,:);
   P(j).orient = AllOrient([-2:0]+j*3,:);
   P(j).weight = AllWeight(j,:);
   % -- modified from original version: JM 30/10/05 --
   P(j).center = AllCenter([-2:0]+j*3,:);
%    P(j).center = [Param.Center]; %one center for both coils
end

[m,n] = size(L);
if(m~=3), % should be 3 x m
   % Old Mosher convention was to give L as m x 3
   % Newer Mathworks convention is for sets of vectors to
   %  be 3 x m (except paradoxically the "patch" command).
   % Error to user, force correction in calling code.
   if Verbose
       bst_message_window('wrap',{'LOCATION GIVEN AS M X 3.',...
               'Please adjust calling code to handle new convention'});
   end

   error('Matrix not given as 3 x n. Correct calling code');
end

G = 0;
for i = 1:length(P),
   G = G + sarvas(L,P(i),Order); % local call below
end

% is there special reference channel considerations?
%  See Channel.mat structure description in the ParameterDescriptions document.
if (refFlag) % Gradient correction is available as well
    
    % read the CTF reference channel information
    meanCenter = mean([Param.Center],2); % mean head center of all of the channels
    % create a temporary parameters file with the same center for all reference channels.
    [RefParam(1:length(irefsens))] = ...
        deal(struct('Center',meanCenter));
    % recursively call
    
    % Forward model on all reference sensors
    % JCM 19-Jun-2002 switched to feval of mfilename
    Gr = feval(mfilename,L,ChannelRef,RefParam,Order,[1:length(irefsens)],[]); % refs now called as MEG
    
    % Apply nth-order gradient correction on good channels only
    
    global ChannelFlag
    if isempty(ChannelFlag)
        ChannelFlag = ones(size(G,1),1); % Take all channels
    end
    
    %Weight by the current nth-order correction coefficients
    try 
        G = G - Channel(1).Comment * Gr; 
    catch
        errordlg(lasterr,...
            'Inconsistency detected in MEG data structure')
        makeuswait('stop')
        return
    end
    
end

clear global Verbose

% ------------- SUBFUNCTIONS, First the simple utilities, then more complicated -----
% -----------------------------------------------------------------------------------
function c = cross(a,b);
% fast and simple, and row major should be faster
% 10/29/99 conversion to row major was slightly faster
%  in the calculation, but overall slower in the transposes needed.
%  retained the column major multiplies below
c = zeros(size(a));
c(1,:) = a(2,:).*b(3,:) - a(3,:).*b(2,:);
c(2,:) = a(3,:).*b(1,:) - a(1,:).*b(3,:);
c(3,:) = a(1,:).*b(2,:) - a(2,:).*b(1,:);





% ----------------------------------------------------------------------------
function k = kronmat(a,b);
% column by column, not element by matrix
k = [a([1 1 1],:) .* b; ...
      a([2 2 2],:) .* b; ...
      a([3 3 3],:) .* b];





% ----------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Local Sarvas functions %%%%%%%%%%%%%%%%%%%%%%%%
% First, the standard interface to above
% If the Order is -1 or 0, then same as before. Order 1 is handled separately

function G = sarvas(L,P,Order);

global Verbose

% Bronzan Sarvas forward model, spherical head
% Order = -1 is current dipole
% Order = 0 is magnetic dipole of BST 2000
% Order = 1 is the new 1st order current multipole
% L is 3 x nL
% 
% P.sensor is 3 x nR,each column a sensor location
% P.orient is 3 x nR, the sensor orientation
% P.center is 3 x nR, the sphere center for each sensor


% January 18, 2002 from sarvas_partial function of 1995

% Used old parameter convention to continue to handle Sylvain's exceptions
%  for the CTF weighting coils

%  if P.center in nonexistant or null, then assumed to be
%  all zeros.

%%%%%%%%%%%%%%%%%%%% which multipolar model to run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch Order
case {-1,0}  %%%%%%%%%%   CURRENT or MAGNETIC DIPOLE %%%%%%%%%%%%%%%%%%%
   
   G = sarvas_dipole(L,P,Order); % BST 2000 function still used here
   
case {1} %%%%%%%%%%%%%%%%%%%%%%%%%  CURRENT MULTIPOLE MODEL %%%%%%%%%%%%%%%%%%555
   % new CME multipole, not MME multipole
   
   % each set of three columns corresponds to a dipole
   % first call the dipolar term
   % SB 02-Jul-2002
   if(size(L,2) > 5*size(P.sensor,2)) & size(P.sensor,2) > 1 & Verbose, %arbitrary definition of a lot
       MESSAGE = 1; % let's show messages
   else
       MESSAGE = 0; % let's keep quiet
   end
   
   if MESSAGE
       bst_message_window('wrap','Initiating the computation of CME modeling . . .')
   end
   
   Gdip = sarvas_dipole(L,P,-1); 
   
   % now call the quadrupole component
   if MESSAGE
       bst_message_window('Refining model . . .') % Not very explicit but gives hint about the progression here - John any better message ?
   end
   
   Gquad = sarvas_partial(L,P); % form the gradient in a local function
   % weights are not applied in this routine, applied separately in the
   %  below forloop
   
   % For example:
   % Gquad(4:6,2) is the gradient of the second coil wrt to the second dipole
   % So we must Vec for each coil to get it into the proper format
   
   nR = size(P.sensor,2); % number of sensors
   nL = size(L,2);  % number of source points
   
   G = zeros(nR,12*nL); % the full gain matrix
   
   for i = 1:nL, % for each dipole
      ndx = [-11:0]+i*12; % indexer for this source
      G(:,ndx(1:3)) = Gdip(:,[-2:0]+i*3); % map next dipole to the first columns
      temp = Gquad(:,[-2:0]+i*3); % next block of partial gradients      
      for j = 1:nR, % for each sensor position, also apply appropriate weight
         G(j,ndx(4:12)) = vec(temp([-2:0]+j*3,:))' * P.weight(j);
      end
   end
   
end

G = G*1e-7; % mu_o over 4 pi







%--------------------------------------------------------------------------
function G = sarvas_dipole(L,P,Order);

global Verbose

% SARVAS MEG Forward Model, spherical head

%  if P.center in nonexistant or null, then assumed to be
%  all zeros.

if(~isfield(P,'center')), % user did not provide
   P.center = []; % initialize to null
end
if(isempty(P.center)), % user gave as null
   P.center = zeros(size(P.sensor));  % set to coordinate origin
end

P.sensor = P.sensor - P.center; % shift sensor coordinates

% SB 03-Sep-2002 
% Now there is the issue of having a sensor array being a mixture of magnetometers and gradiometers.
% Magnetometers are referred as pseudo-gradiometers: the corresponding .Loc field of the Channel array of structures
% is still 3x2. For magnetometers though, the 2nd colum is filled with zeros (i.e.: [0 0 0]').
% Same story holds for the .Orient field.
% The full gain matrices are computed for these channels located at [0 0 0] but when the orientation vector is applied
% the net field is set to 0. Therefore, the second call to sarvas_dipole produces a null field which is substracted frim the field 
% from the magnetometer (i.e the fisrt coil of the pseudo-gradiometer).
% Calculations fail (divide by zero) when head center is also at [0 0 0]. 
% I'm therefore testing this out here and fix things by virtually translating the virtual coils located at [0 0 0]
% to [1 1 1] (arbitrary). The net field will still be null when the orientation is applied anyway.
% Any more elegant fix is welcome at this point.
iMag = find(norcol(P.sensor)==0); % Indices of channels located at P.center.
if ~isempty(iMag)
    P.sensor(:,iMag) = repmat([1 1 1]',1,length(iMag)); % Move them away (arbitrary location).
end

nR = size(P.sensor,2); % number of sensors
nL = size(L,2);  % number of source points

Rn2 = sum(P.sensor.^2,1); % distance to sensor squared
Rn = sqrt(Rn2); % distance

if (nR >= nL), % more sensors than dipoles
   if(Order == 1),
      G = zeros(nR,12*nL); % gain matrix
   else    
      G = zeros(nR,3*nL);  % gain matrix
   end
   
   for Li = 1:nL,
      Lmat = L(:,Li+zeros(1,nR)); % matrix of location repeated
      Lmat = Lmat - P.center; % each center shifted relative to its center
      D = P.sensor - Lmat;  % distance from souce to sensors
      Dn2 = sum(D.^2,1); % distance squared
      Dn = sqrt(Dn2);  % distance
      R_dot_D = sum(P.sensor .* D);  % dot product of sensor and distance
      R_dot_Dhat = R_dot_D ./ Dn;  % dot product of sensor and distance
      
      F = Dn2 .* Rn + Dn .* R_dot_D;  % Sarvas' function F
      
      GF_dot_o = Dn2 .* sum(P.sensor.*P.orient) ./ Rn + ...
         (2 * Rn + R_dot_Dhat) .* sum(D.*P.orient) + ...
         Dn .* sum((D+P.sensor).*P.orient);
      
      tempF = GF_dot_o ./ F.^2;
      
      if(Order == -1), % current dipole model
         temp = cross(Lmat,P.orient) ./ F([1 1 1],:) - ...
            cross(Lmat,P.sensor) .* tempF([1 1 1],:);
         G(:,Li*3+[-2 -1 0]) = temp';
      elseif(Order == 0) % magnetic dipole model
         temp = P.sensor .* tempF([1 1 1],:) - P.orient ./ F([1 1 1],:);
         G(:,Li*3+[-2 -1 0]) = temp';
      elseif(Order == 1),  % 1st order multipole
         % first the dipole
         temp_m = P.sensor .* tempF([1 1 1],:) - P.orient ./ F([1 1 1],:);
         % then the quadrupole
         temp1 = -(2*Rn + R_dot_Dhat + Dn);
         temp2 = -sum(D.*P.orient) ./ Dn;
         temp3 = -(2*sum(P.sensor.*P.orient)./Rn + sum((D+P.sensor).*P.orient)./Dn - ...
            sum(D.*P.orient).*R_dot_D./(Dn2.*Dn));
         
         GGpF_dot_o = temp1([1 1 1],:) .* P.orient + ...
            temp2([1 1 1],:).*P.sensor + temp3([1 1 1],:) .* D;
         temp1 = -(2*Rn + R_dot_Dhat);
         GpF = temp1([1 1 1],:) .* D - Dn([1 1 1],:) .* P.sensor;
         temp1 = 1 ./ F.^2;
         temp2 = 2*GF_dot_o./F;
         temp_q = temp1(ones(1,9),:) .* (kronmat(GGpF_dot_o,P.sensor) + ...
            kronmat(GpF,P.orient - temp2([1 1 1],:).*P.sensor));
         G(:,Li*12+[-11:0]) = [temp_m;temp_q]'; 
      end
      
   end
   
else  % more dipoles than sensors nL > nR
   if(Order == 1)
      G = zeros(12*nL,nR);  % 1st order multipole gain matrix transposed
   else
      G = zeros(3*nL,nR);  % gain matrix transposed
   end
   
   % if there are a lot of dipoles, let's watch on the screen
   % JCM 18-Jun-2002 no more waitbar, use message window
   if (nL > 5*nR) & nR > 1 & Verbose, %arbitrary definition of a lot
      MESSAGE = 1; % let's show messages
      %       bst_message_window('wrap',sprintf(...
      %          'Making order %.0f matrix of %.0f sensors x %.0f sources',Order,nR,nL));
      %       bst_message_window('append','Making first sensor . . .'); % prepare for overwrite
      bst_message_window('overwrite',sprintf(...
          'Making Current Dipole matrix of %.0f sensors x %.0f sources',nR,nL));
      bst_message_window('append','Making first sensor . . .'); % prepare for overwrite   else
  else % SB 03-Sep-2002 : else was probably missing 
      MESSAGE = Verbose; % let's be quiet
  end
  
   for Ri = 1:nR,
      if(MESSAGE), % want to show the user progress?
         if(~rem(Ri,30)), % every tenth sensor
            bst_message_window('overwrite',sprintf('Progress report:....... %.0f of %.0f . . .',Ri,nR));
         end
      end
      Rmat = P.sensor(:,Ri+zeros(1,nL)); % matrix of sensor repeated
      Omat = P.orient(:,Ri+zeros(1,nL)); % orientations
      Lmat = L - P.center(:,Ri+zeros(1,nL)); % shift centers to this coordinate
      
      D = Rmat - Lmat;
      Dn2 = sum(D.^2,1); % distance squared
      Dn = sqrt(Dn2);  % distance
      R_dot_D = sum(Rmat .* D);  % dot product of sensor and distance
      R_dot_Dhat = R_dot_D ./ Dn;  % dot product of sensor and distance
      
      F = Dn2 * Rn(Ri) + Dn .* R_dot_D;  % Sarvas' function F
      
      GF_dot_o = Dn2 * sum(P.sensor(:,Ri).*P.orient(:,Ri)) / Rn(Ri) + ...
         (2 * Rn(Ri) + R_dot_D ./ Dn) .* sum(D.*Omat) + ...
         Dn .* sum((D+Rmat).*Omat);
      
      tempF = GF_dot_o ./ F.^2;
      
      if(Order == -1), % current dipole model
         temp = cross(Lmat,Omat) ./ F([1 1 1],:) - ...
            cross(Lmat,Rmat) .* tempF([1 1 1],:);
      elseif(Order == 0) % magnetic dipole model
         temp = Rmat .* tempF([1 1 1],:) - Omat ./ F([1 1 1],:);
      elseif(Order == 1),  % 1st order multipole
         % first the dipole
         temp_m = Rmat .* tempF([1 1 1],:) - Omat ./ F([1 1 1],:);
         % then the quadrupole
         temp1 = -(2*Rn(Ri) + R_dot_Dhat + Dn);
         temp2 = -sum(D.*Omat) ./ Dn;
         temp3 = -(2*sum(P.sensor(:,Ri).*P.orient(:,Ri))./Rn(Ri) + sum((D+Rmat).*Omat)./Dn - ...
            sum(D.*Omat).*R_dot_D./(Dn2.*Dn));
         
         GGpF_dot_o = temp1([1 1 1],:) .* Omat + ...
            temp2([1 1 1],:).*Rmat + temp3([1 1 1],:) .* D;
         temp1 = -(2*Rn(Ri) + R_dot_Dhat);
         GpF = temp1([1 1 1],:) .* D - Dn([1 1 1],:) .* Rmat;
         temp1 = 1 ./ F.^2;
         temp2 = 2*GF_dot_o./F;
         temp_q = temp1(ones(1,9),:) .* (kronmat(GGpF_dot_o,Rmat) + ...
            kronmat(GpF,Omat - temp2([1 1 1],:).*Rmat));
         temp = [temp_m;temp_q];
      else % unimplemented order
         disp('SARVAS: Unimplemented source order');
         G = [];
         return
      end    
      
      G(:,Ri) = temp(:);
   end
   
   G = G';
   
end

if(isfield(P,'weight')),
   Weights = P.weight(:); %make sure column
   % scale each row by its appropriate weight
   G = Weights(:,ones(1,size(G,2))) .* G;
end





% ----------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%   Partial calculation for the sarvas %%%%%%%%%%%%%%%%%%%

function D = sarvas_partial(L,P);
%SARVAS_PARTIAL Calculate the partial of the Sarvas Formula w.r.t L
% function D = sarvas_partial(L,P);
% For dipole locations in L and sensor information in R, calculate the
% partial of the Sarvas formula (dipole in a sphere, arbitrary sensor
% orientation) with respect to the dipole location.  The result is a matrix
% D of partials information.  
% If there are M sensors and P dipoles, then D
% is 3*M by 3*P.  Let the corresponding moments Q be 3 by M.
% Then partial = D * blk_diag(Q,1) is 3*M by P.  Each column of partial
% corresponds to a different dipole in L.  Each set of three rows
% corresponds to a different sensor location.  Thus partial(4:6,2) is the
% partial of the Sarvas formula at the second sensor location, with respect
% to the second dipole.
%
% Used to calculate the Cramer-Rao Lower Bounds
% structure P has fields .weight, length of the number of sensors, which
%  is applied to each row of the gain matrix. Weight is usually 1 or -1
%  Field .center is 3 x number of sensors, head center for each coil.
%  P.sensor and P.orient are the location and orientation information, resp.
%  P.weight is not applied, must be applied separately, see above calling code.

% Author: John C. Mosher, Ph.D.
% Los Alamos National Laboratory
% Los Alamos, NM 87545
% email: mosher@LANL.Gov

% August 15, 1995 author
% 20 Feb 2002 JCM extensive I/O changes to adapt into os_meg, handle arbitrary 
%                 head center

L = mby3check(L,0);		% matrices are now 3 x <>, but don't warn if 3 x 3

numSens = size(P.sensor,2); % number of sensors
numDip = size(L,2);  % number of source points

% use some older notation for programming reuse below
R = [P.sensor]; % sensor locations
Rs = [P.orient]; % sensor orientations

if(~isfield(P,'center')), % user did not provide
   P.center = []; % initialize to null
end
if(isempty(P.center)), % user gave as null
   P.center = zeros(size(R));  % set to coordinate origin
end

% each channel coil is shifted to its relative location
R = R - P.center; % shift sensor coordinates

% SB 03-Sep-2002 
% Now there is the issue of having a sensor array being a mixture of magnetometers and gradiometers.
% Magnetometers are referred as pseudo-gradiometers: the corresponding .Loc field of the Channel array of structures
% is still 3x2. For magnetometers though, the 2nd colum is filled with zeros (i.e.: [0 0 0]').
% Same story holds for the .Orient field.
% The full gain matrices are computed for these channels located at [0 0 0] but when the orientation vector is applied
% the net field is set to 0. Therefore, the second call to sarvas_dipole produces a null field which is substracted frim the field 
% from the magnetometer (i.e the fisrt coil of the pseudo-gradiometer).
% Calculations fail (divide by zero) when head center is also at [0 0 0]. 
% I'm therefore testing this out here and fix things by virtually translating the virtual coils located at [0 0 0]
% to [1 1 1] (arbitrary). The net field will still be null when the orientation is applied anyway.
% Any more elegant fix is welcome at this point.
iMag = find(norcol(R)==0); % Indices of channels located at P.center.
if ~isempty(iMag)
    R(:,iMag) = repmat([1 1 1]',1,length(iMag)); % Move them away (arbitrary location).
end

Rn = colnorm(R); 		% distance to sensor 
iRn = 1../Rn;			% inverse distance
o3 = ones(3,1);			% col vector of three ones

% three by three matrices per sensor and dipole
D = zeros(numSens*3,3*numDip);

if(size(L,2) > 5*size(R,2)) & size(R,2) > 1 & Verbose, %arbitrary definition of a lot
    MESSAGE = Verbose; % let's show messages
    bst_message_window({...
            'Computing quadrupolar moments',...
            sprintf('for sources %.0f of %.0f . . .',100,numDip)...
        });
    
else
    MESSAGE = Verbose; % let's be quiet
end

for Dip = 1:numDip		% foreach dipole,

    %################ main loop ########################
  
   % if there are a lot of dipoles, let's watch on the screen
   % SB 02-Jul-2002 no more waitbar, use message window
   
   if(MESSAGE), % want to show the user progress?
       if(~rem(Dip,100)), % every tenth source
           bst_message_window('overwrite',sprintf('for sources %.0f of %.0f . . .',Dip,numDip));
       end
   end
   
   ThisDipole = L(:,Dip);		% next dipole
   Lmat = repmat(ThisDipole,1,numSens); % repeat this dipole for all sensors
   Lmat = Lmat - P.center; % shift each location relative to the sensor's center
   
   % let "a" be the same as Sarvas' "a", a = sensor - dipole.
   a = R - Lmat;
   an = colnorm(a);		% norm of each a
   ian=1../an;			% inverse of norm a
   ian3 = ian.^3;			% inverse cubed of norm a
   
   aDotR = sum(a.*R,1);
   
   % Form F
   
   F = (an .* Rn + aDotR) .* an;
   iF = 1../F;
   
   % From gradient F
   
   tmp1 = an.^2 .* iRn + aDotR .* ian + 2*an + 2*Rn;
   tmp1 = tmp1(o3,:).* R;
   tmp2 = an + 2*Rn + aDotR .* ian;
   tmp2 = tmp2(o3,:) .* Lmat;
   
   gradF = tmp1 - tmp2;
   
   % take partials of (Rs dot grad(F)) wrt r_q.  Result is 3 by 1 per sensor
   
   tmp1 = aDotR.*ian3 - 2*(iRn+ian);
   tmp1 = (tmp1(o3,:) .* a) - (R .* ian(o3,:));
   tmp1 = tmp1 .* (o3*sum(Rs.*R,1));
   
   tmp2 = an + 2*Rn + ian.*aDotR;
   tmp2 = tmp2(o3,:) .* Rs;
   
   tmp3 = aDotR.*ian3 - ian;
   tmp3 = tmp3(o3,:).*a - (R .* ian(o3,:));
   tmp3 = tmp3 .* (o3*sum(Rs.*Lmat,1));
   
   partRsGradF = tmp1 - tmp2 - tmp3;
   
   
   % take partial of F wrt r_q
   
   tmp1 = aDotR .* ian;
   tmp1 = tmp1(o3,:).*a;
   
   partF = -2*Rn(o3,:).*a - an(o3,:).*R - tmp1;
   
   
   % now take partial of ((gradF dot s) / F^2).
   
   tmp1 = F(o3,:) .* partRsGradF;
   
   tmp2 = 2*sum(gradF.*Rs,1);
   tmp2 = tmp2(o3,:) .* partF;
   
   partGradFdotRsOverF2 = (tmp1 - tmp2) .* iF(o3,:).^3;
   
   
   % now take partial of inverse of F
   
   partInvF = -partF .* iF(o3,:).^2;
   
   
   %  Now we are ready to generate the partial of the Sarvas' formula wrt dipole
   %  location.  The result is a 3 by 1 per sensor location; however, we want
   %  to separate out the dipole moment q.  So our "partials matrix" D is 3 x 3
   %  per sensor location.  We will concatenate into a 3*M by 3 matrix for each
   %  dipole.
   
   % Rs cross q divide by F
   % Want the cross product matrix tensor
   
   tmp1 = blk_diag(cross_mat(Rs),3); % each set of three columns is a tensor
   tmp1 = tmp1 .* kron(iF,ones(3)); % divide each tensor by F
   tmp1 = blk_lex(tmp1,3); 	% now each set of three rows is a tensor
   
   % L cross Rs dot q time partInvF
   % Want the direct (outer) product of partInvF and (L cross Rs)
   
   tmp2 = cross(Lmat,Rs); % cross products
   % each set of three columns is the outer product
   tmp2 = partInvF * blk_diag(tmp2,1)';
   tmp2 = blk_lex(tmp2,3);		% now each set of three rows is a tensor
   
   % (R cross q)*(gradF dot Rs)/F^2
   % Want scalar times the cross product tensor
   
   tmp3 = sum(gradF .* Rs,1) .* iF.^2;
   tmp3 = kron(tmp3,ones(3));	% repeat scalar for each submatrix
   tmp3 = tmp3 .* blk_diag(cross_mat(R),3);
   tmp3 = blk_lex(tmp3,3);		% vertically stacked now
   
   % ((L cross R) dot q) * partial(gradF dot Rs over F^2)
   % Want direct (outer) product of partGradFdotRsOverF2 and (L cross R)
   
   tmp4 = cross(Lmat,R);
   tmp4 = partGradFdotRsOverF2 * blk_diag(tmp4,1)';
   tmp4 = blk_lex(tmp4,3);
   
   % Now combine into the appropriate columns of D.
   
   D(:,(Dip-1)*3 + [1:3]) = tmp1 + tmp2 - tmp3 - tmp4;
   
   %################ end main loop ########################
end 				% next dipole
