function [Channel,B] = spm_eeg_rdata_FIF_channels(Fdata);

% Function that returns channel info for various types of Neuromag coils
% in BrainStorm format (made separate from spm_eeg_rdata_FIF in case 
% such channel info is needed for Brainstorm directly, without using SPM)
%
% Based on "importdata.m" function of Brainstorm, but updated with Appendix B
% of Elekta Neuromag Source Modelling Software Usersguide v5.4 (Xfit-5.4.pdf)
%
% NOTES: 
%    - Requires the Fiff Access toolbox written by Kimmo Uutela
%      http://www.kolumbus.fi/~w132276/programs/meg-pd/
%    - Fiff Access does not handle long int, so may need to convert 
%      to short int or float first (eg using MaxFilter)
%
% Rik Henson (5/06/07)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This cell array of structures contain the information on every type of 
% coil that is used in the Neuromag/BTi systems

Coil = cell(4001,1);

Coil{2}.Description = 'Neuromag-122 planar gradiometer';
Coil{2}.n = 2; % Number of integration points
Coil{2}.xyz = [8.1 0 0 ; -8.1 0 0]'; % Locations of Integration points for the magnetic field (mm)
Coil{2}.weights = [1/16.2 -1/16.2]; % weights to be applied at each integration point  

Coil{3012}.Description = 'VectorView type 1 planar gradiometer';
Coil{3012}.n = 2; % Number of integration points
Coil{3012}.xyz = [8.4 0 0.3 ; -8.4 0 0.3]'; % Locations of Integration points for the magnetic field (mm)
Coil{3012}.weights = [1/16.8 -1/16.8]; % weights to be applied at each integration point  

Coil{3013}.Description = 'VectorView type 2 planar gradiometer';
Coil{3013}.n = 2; % Number of integration points
Coil{3013}.xyz = [8.4 0 0.3 ; -8.4 0 0.3]'; % Locations of Integration points for the magnetic field (mm)
Coil{3013}.weights = [1/16.8 -1/16.8]; % weights to be applied at each integration point  

Coil{3014}.Description = 'VectorView type 3 planar gradiometer';
Coil{3014}.n = 2; % Number of integration points
Coil{3014}.xyz = [8.4 0 0.3 ; -8.4 0 0.3]'; % Locations of Integration points for the magnetic field (mm)
Coil{3014}.weights = [1/16.8 -1/16.8]; % weights to be applied at each integration point  

Coil{3022}.Description = 'VectorView type 1 magnetometer';
Coil{3022}.n = 4; % Number of integration points
%Coil{3022}.xyz = [12.9 12.9 0.3; 12.9 -12.9 0.3; -12.9 12.9 0.3; -12.9 -12.9 0.3]'; % Locations of Integration points for the magnetic field (mm)
Coil{3022}.xyz = [6.45 6.45 0.3; 6.45 -6.45 0.3; -6.45 6.45 0.3; -6.45 -6.45 0.3]'; % Locations of Integration points for the magnetic field (mm)
Coil{3022}.weights = [1/4 1/4 1/4 1/4]; % weights to be applied at each integration point  

Coil{3023}.Description = 'VectorView type 2 magnetometer';
Coil{3023}.n = 4; % Number of integration points
%Coil{3023}.xyz = [12.9 12.9 0.3 ;12.9 -12.9 0.3; -12.9 12.9 0.3; -12.9 -12.9 0.3]'; % Locations of Integration points for the magnetic field (mm)
Coil{3023}.xyz = [6.45 6.45 0.3; 6.45 -6.45 0.3; -6.45 6.45 0.3; -6.45 -6.45 0.3]'; % Locations of Integration points for the magnetic field (mm)
Coil{3023}.weights = [1/4 1/4 1/4 1/4]; % weights to be applied at each integration point  

Coil{3024}.Description = 'VectorView type 3 magnetometer';
Coil{3024}.n = 4; % Number of integration points
Coil{3024}.xyz = [5.25 5.25 0.3; 5.25 -5.25 0.3; -5.25 5.25 0.3; -5.25 -5.25 0.3]'; % Locations of Integration points for the magnetic field (mm)
Coil{3024}.weights = [1/4 1/4 1/4 1/4]; % weights to be applied at each integration point  

Coil{2000}.Description = 'Ideal point magnetometer';
Coil{2000}.n = 1; % Number of integration points
Coil{2000}.xyz = [0 0 0];
Coil{2000}.weights = [1]; % weights to be applied at each integration point  

Coil{4001}.Description = 'Magnes WH magnetometer';
Coil{4001}.n = 4; % Number of integration points
Coil{4001}.xyz = [5.75 5.75 0; -5.75 5.75 0;...
-5.75 -5.75 0; 5.75 -5.75 0]';
Coil{4001}.weights = [1/4 1/4 1/4 1/4]; % weights to be applied at each integration point  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ChanName,ChanType,ChanNum,CoilType,T] = channames(Fdata);

MegChan = find(ChanType==1);

nchan = length(MegChan);
 
[Channel(1:nchan)] = deal(struct('Loc',[],'Orient',[],'Comment',[],'Weight',[],'Type',[],'Name','','SensType',''));
    
for c = 1:nchan

    chan = MegChan(c);

% This part from Brainstorm importdata

    Channel(chan).CoilLoc = T{chan}*[Coil{CoilType(chan)}.xyz/1000; ones(1,Coil{CoilType(chan)}.n)];
    Channel(chan).CoilLoc = Channel(chan).CoilLoc(1:3,:);
    Channel(chan).CoilOrient = T{chan}*[[0 0 1]' * ones(1,Coil{CoilType(chan)}.n); zeros(1,Coil{CoilType(chan)}.n)];
    Channel(chan).CoilOrient = Channel(chan).CoilOrient(1:3,:);
    Channel(chan).Comment = ChanName(chan,:);
    Channel(chan).CoilWeight = Coil{CoilType(chan)}.weights;  
    Channel(chan).Type = 'MEG';
    Channel(chan).Name = strrep(Channel(chan).Comment,'MEG ','');  

% Next is how Brainstorm (bst_headmodeller) wants the data passed (in OPTIONS)...

    if size(Channel(chan).CoilLoc,2) == 2
% Hack for planar gradiometers

	Channel(chan).SensType = 'Gradiometer';
	Channel(chan).Loc = [Channel(chan).CoilLoc(:,1); Channel(chan).CoilLoc(:,2)];
	Channel(chan).Weight = Channel(chan).CoilWeight;

	% This is orientation of COIL (ie loop)
        Channel(chan).Orient = [Channel(chan).CoilOrient(:,1); Channel(chan).CoilOrient(:,2)];		

	% This is orientation of (planar) GRADIENT (ie difference measured by coils, which is X-direction)
	Channel(chan).GradOrient = T{chan}*[[1 0 0]' * ones(1,Coil{CoilType(chan)}.n); zeros(1,Coil{CoilType(chan)}.n)];
	Channel(chan).GradOrient = Channel(chan).GradOrient(1:3,1);		

    elseif size(Channel(chan).CoilLoc,2) == 4
% Hack for magnetometers
	Channel(chan).SensType = 'Magnetometer';

% Comments in spm_bst_headmodeler require NaN for "second coil" of 
% magnetometers (though spm_eeg_inv_BSTfwdsol actually handles this)
	Channel(chan).Loc = [mean(Channel(chan).CoilLoc,2); ones(3,1)*NaN];
	Channel(chan).Orient = [Channel(chan).CoilOrient(1:3,1); ones(3,1)*NaN];
	Channel(chan).Weight = [sum(Channel(chan).CoilWeight) 0];
	Channel(chan).GradOrient = Channel(chan).Orient(1:3);		

    else
	error(sprintf('Sensor %d: Unknown coil type',chan))
    end
end

B.chans = ChanName;
B.chtypes = ChanType;
B.chnums = ChanNum;
B.coiltypes = CoilType;
B.chanfilt = MegChan;
B.T = T;
