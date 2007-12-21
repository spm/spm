function [fid, sens] = spm_eeg_inv_ReadPolhemus(Fname_pol,skip,figflag);

%=======================================================================
% Reads Polhemus files:
%   either sensor file or headshape file or both
%
% FORMAT [fid, sens] = spm_eeg_inv_ReadPolhemus(Fname.pol, skip, figflag)
% Input:
% Fname_pol - Polhemus ASCII file containing sensor locations (cm)
%             (headshape can also be considered here instead of sensors)
% skip      - channels to skip
% figflag   - display the point locations (1) or not (0) (default: 0)
%
% Output:
% fid       - fiducial         locations (mm) in rows
% sens      - sensor/headshape locations (mm) in rows
%
% IMPORTANT: Note that Polhemus data files should be -ASCII files with
% extension .pol
% It is assumed that the .pol file contains the location (cm) of fiducials
% (sampled twice), followed by the location of the sensors.  In some
% instances the first few channel locations may pertain to reference
% channels; the skip variable allows these to be skipped if necessary.
% The fiducial locations are flaged with the strings 'NZ','LE' and 'RE'; 
% indicating the Nasion, left and right eare respectively
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_ReadPolhemus.m 1039 2007-12-21 20:20:38Z karl $



% checks and assigments
%--------------------------------------------------------------------------
try, skip;    catch, skip    = 0; end
try, figflag; catch, figflag = 0; end

[pth,nam,ext] = fileparts(Fname_pol);
if ~strcmp(ext,'.pol')
    warndlg(sprintf('Wrong input file format\n'));
    return
end

if figflag
    F = spm_figure('GetWin','Graphics');
    figure(F);  clf
end


% --- READ Polhemus Sensor + fiducial locations ---
%==========================================================================
try
    file = textread(Fname_pol,'%s');
catch
    file = textread(fullfile(pwd,[nam ext]),'%s');
end
% remove zeros at the end
bool = 0;
while bool == 0
    if strcmp(file{end},'0')
        file = file(1:end-1);
    else
        bool = 1;
    end
end

% read fiducials
%--------------------------------------------------------------------------
NZ   = [];
LE   = [];
RE   = [];
temp = 0;
nl   = 1;

while temp == 0
    if strcmp(file{nl},'NZ')
        NZ = [NZ ; str2num(file{nl+1}) str2num(file{nl+2}) str2num(file{nl+3}) ];
        nl = nl + 4;
    elseif strcmp(file{nl},'LE') | strcmp(file{nl},'OG')
        LE = [LE ; str2num(file{nl+1}) str2num(file{nl+2}) str2num(file{nl+3}) ];
        nl = nl + 4;
    elseif strcmp(file{nl},'RE') | strcmp(file{nl},'OD')
        RE = [RE ; str2num(file{nl+1}) str2num(file{nl+2}) str2num(file{nl+3}) ];
        nl = nl + 4;
    else
        temp = 1;
    end
end

% convert from cm to mm
%--------------------------------------------------------------------------
NZ    = mean(NZ,1); LE = mean(LE,1); RE = mean(RE,1);
fid   = 10*[NZ; LE; RE]; 

% read sensor lcoations or headshape locations
%--------------------------------------------------------------------------
sens  = [];
start = nl + skip*3;
for i = start:3:length(file)
    sens = [sens; str2num(file{i}) str2num(file{i+1}) str2num(file{i+2})];
end

% convert from cm to mm
%--------------------------------------------------------------------------
sens   = 10*sens;


% Display
%--------------------------------------------------------------------------
if figflag

    h_sens = plot3(sens(:,1),sens(:,2),sens(:,3),'or');
    set(h_sens,'MarkerFaceColor','r','MarkerSize',12,'MarkerEdgeColor','k');
    hold on
    h_fid1 = plot3(fid(:,1),fid(:,2),fid(:,3),'om');
    set(h_fid1,'MarkerFaceColor','m','MarkerSize',12,'MarkerEdgeColor','k');

    % cameramenu
    axis off;
    axis image;
    rotate3d
    drawnow
end
