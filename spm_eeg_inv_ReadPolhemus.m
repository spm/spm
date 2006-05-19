function [varargout] = spm_eeg_inv_ReadPolhemus(varargin);

%=======================================================================
% Reads Polhemus files (for EEG mainly):
%   either sensor file or headshape file or both
%   if .pol files are not provided as inputs, user will be prompted
%   .mat files required by SPM for data registration will be created
%   with coordinates in mm
%
% FORMAT D = spm_eeg_inv_ReadPolhemus(S)
% Input:
% S		    - input data struct (optional)
% Output:
% D         - same data structure containing the created file names for
%           data registration
%
% FORMAT [Fname_fid.mat , Fname_sens.mat] = spm_eeg_inv_ReadPolhemus(Fname.pol, figflag)
% Input:
% Fname_pol - Polhemus ASCII file containing sensor locations
%           (headshape can also be considered here instead of sensors)
% figflag   - display the point locations (1) or not (0) (default: 0)
% Output:
% Fname...mat - input Polhemus data have been save in matrix (Matlab) format
%
% FORMAT [Fname_fid.mat , Fname_sens.mat , Fname_fid2.mat , Fname_hsp.mat] = spm_eeg_inv_ReadPolhemus(Fname_sens.pol , Fname_hsp.pol, figflag)
% Input:
% Fname_sens.pol - Polhemus ASCII file containing sensor (& fiducials)
% Fname_hsp.pol  - Polhemus ASCII file containing headshape (& fiducials)
% figflag   - display the point locations (1) or not (0) (default: 0)
%
% Output:
% Fname....mat
%
% IMPORTANT: Note that Polhemus data files should be -ASCII files with
% extension .pol
%=======================================================================
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Jeremie Mattout
% $Id: spm_eeg_inv_ReadPolhemus.m 539 2006-05-19 17:59:30Z Darren $

spm_defaults

Filenames = {};
def_figflag = 0;

if nargout == 0 | nargout == 1

    try
        D = varargin{1};
    catch
        D = spm_select(1, '.mat', 'Select EEG/MEG mat file');
        D = spm_eeg_ldata(D);
    end

    if ~isfield(D,'inv')
        error(sprintf('No inverse structure has been created for this data set\n'));
    end
    
    val = length(D.inv);
   
    try
        Isens = D.channels.eeg;
    catch
        Isens = [];
    end
    
    if nargin > 1
        figflag = varargin{2};
    else
        figflag = def_figflag;
    end

elseif nargout == 2
    
    if nargin >= 1 & nargin <= 2
        Filenames{1} = varargin{1};
        [pth,nam,ext] = fileparts(Filenames{1});
        if ~strcmp(ext,'.pol')
            error(sprintf('Wrong input file format\n'));
        end
        if nargin == 2
            figflag = varargin{2};
        else
            figflag = def_figflag;
        end
    else
        error(sprintf('Wrong input arguments: one .pol files required\n'));    
    end
    
elseif nargout == 4
    
    if nargin >= 2 & nargin <= 3
        Filenames{1} = varargin{1};
        [pth,nam,ext] = fileparts(Filenames{1});
        if ~strcmp(ext,'.pol')
            error(sprintf('Wrong input file format: .pol required\n'));
        end
        Filenames{2} = varargin{2};
        [pth,nam,ext] = fileparts(Filenames{2});
        if ~strcmp(ext,'.pol')
            error(sprintf('Wrong input file format: .pol required\n'));
        end
        figflag = def_figflag;
        if nargin == 3
            figflag = varargin{3};
        else
            figflag = def_figflag;
        end
    else
        error(sprintf('Wrong input arguments: two .pol files required\n'));    
    end
         
elseif nargout > 4 | nargin > 3
        
    error(sprintf('Wrong output arguments\n'));
    
end

if isempty(intersect(figflag,[0 1]))
    error(sprintf('Display flag should be either 0 (no display) or 1 (display)'));
end

switch length(Filenames)
    case 0
        Filenames{1} = spm_select(1, '.pol', 'Select sensor Polhemus file');
        Filenames{2} = spm_select(1, '.pol', 'Select headshape Polhemus file');
    case 1
        HSPflag = spm_input('Download Head Shape?','+1','Yes|No',[1 2]);
        if HSPflag == 1
            Filenames{2} = spm_select(1, '.pol', 'Select headshape Polhemus file');
        else
            Filenames{2} = '';
        end
end

if figflag
    F = findobj('Tag', 'Graphics');
    
    if isempty(F)
        F = spm_figure;
    end
    
    figure(F);
    clf
end


% --- READ Polhemus Sensor + fiducial locations ---
file = textread(Filenames{1},'%s');
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
NZ = [];
LE = [];
RE = [];
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
NZ = mean(NZ); LE = mean(LE); RE = mean(RE);
fid = 10*[NZ ; LE ; RE]; % convert from cm to mm
try
    [pth,nam,ext] = fileparts(D.inv{val}.mesh.sMRI);
    fidname = [nam '_fid_EEG_sens.mat'];
    D.inv{val}.datareg.fid = fidname;
    if spm_matlab_version_chk('7') >=0
        save(fullfile(D.path, fidname), '-V6', 'fid');
    else
    	save(fullfile(D.path, fidname), 'fid');
    end
catch
    fidname = 'fid_EEG_sens.mat';
    if spm_matlab_version_chk('7') >=0
        save(fidname, '-V6', 'fid');
    else
    	save(fidname, 'fid');
    end
end

sens = [];
for i = nl:3:length(file)
    sens = [sens ; str2num(file{i}) str2num(file{i+1}) str2num(file{i+2})];
end
sens = 10*sens; % convert from cm to mm
try
    % select channels
    if ~isempty(Isens) & length(Isens) <= length(sens)
        sens = sens(end-length(Isens)+1:end,:);
    end    
    [pth,nam,ext] = fileparts(D.inv{val}.mesh.sMRI);
    sensname = [nam '_sensors_EEG_sens.mat'];
    D.inv{val}.datareg.fid = sensname;
    if spm_matlab_version_chk('7') >=0
        save(fullfile(D.path, sensname), '-V6', 'sens');
    else
    	save(fullfile(D.path, sensname), 'sens');
    end
catch
    sensname = 'sensors_EEG_sens.mat';
    if spm_matlab_version_chk('7') >=0
        save(sensname, '-V6', 'sens');
    else
    	save(sensname, 'sens');
    end
end
clear file


% --- READ Polhemus HeadShape + fiducial locations ---
if ~isempty(Filenames{2})
    file = textread(Filenames{2},'%s');
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
    NZ2 = [];
    LE2 = [];
    RE2 = [];
    temp = 0;
    nl   = 1;

    while temp == 0
        if strcmp(file{nl},'NZ')
            NZ2 = [NZ2 ; str2num(file{nl+1}) str2num(file{nl+2}) str2num(file{nl+3}) ];
            nl = nl + 4;
        elseif strcmp(file{nl},'LE') | strcmp(file{nl},'OG')
            LE2 = [LE2 ; str2num(file{nl+1}) str2num(file{nl+2}) str2num(file{nl+3}) ];
            nl = nl + 4;
        elseif strcmp(file{nl},'RE') | strcmp(file{nl},'OD')
            RE2 = [RE2 ; str2num(file{nl+1}) str2num(file{nl+2}) str2num(file{nl+3}) ];
            nl = nl + 4;
        else
            temp = 1;
        end
    end
    NZ2 = mean(NZ2); LE2 = mean(LE2); RE2 = mean(RE2);
    fid2 = 10*[NZ2 ; LE2 ; RE2]; % convert from cm to mm
    try
        [pth,nam,ext] = fileparts(D.inv{val}.mesh.sMRI);
        fidname2 = [nam '_fid_EEG_hsp.mat'];
        if spm_matlab_version_chk('7') >=0
            save(fullfile(D.path, fidname2), '-V6', 'fid2');
        else
        	save(fullfile(D.path, fidname2), 'fid2');
        end
    catch
        fidname2 = 'fid_EEG_hsp.mat';
        if spm_matlab_version_chk('7') >=0
            save(fidname2, '-V6', 'fid2');
        else
        	save(fidname2, 'fid2');
        end
    end
    
    hsp = [];
    for i = nl:3:length(file)
        hsp = [hsp ; str2num(file{i}) str2num(file{i+1}) str2num(file{i+2})];
    end
    hsp = 10*hsp;
    try
        [pth,nam,ext] = fileparts(D.inv{val}.mesh.sMRI);
        hspname = [nam '_headshape_EEG.mat'];
        D.inv{val}.datareg.fid = hspname;
        if spm_matlab_version_chk('7') >=0
            save(fullfile(D.path, hspname), '-V6', 'hsp');
        else
        	save(fullfile(D.path, hspname), 'hsp');
        end
    catch
        hspname = 'headshape_EEG.mat';
        if spm_matlab_version_chk('7') >=0
            save(hspname, '-V6', 'hsp');
        else
        	save(hspname, 'hsp');
        end
    end
end


% Save results
if nargout == 1
    if spm_matlab_version_chk('7') >= 0
        save(fullfile(D.path, D.fname), '-V6', 'D');
    else
     	save(fullfile(D.path, D.fname), 'D');
    end
    varargout{1} = D;
elseif nargout == 2
    varargout{1} = fidname;
    varargout{2} = sensname;
elseif nargout == 4
    varargout{1} = fidname;
    varargout{2} = sensname;
    varargout{3} = fidname2;
    varargout{4} = hspname;
end
    
        
 
% Display
if figflag
    
    % Display sensors
    h_sens = plot3(sens(:,1),sens(:,2),sens(:,3),'or');
    set(h_sens,'MarkerFaceColor','r','MarkerSize',12,'MarkerEdgeColor','k');
    hold on
    h_fid1 = plot3(fid(:,1),fid(:,2),fid(:,3),'om');
    set(h_fid1,'MarkerFaceColor','m','MarkerSize',12,'MarkerEdgeColor','k');
    
    % Display head-shape if any
    if ~isempty(Filenames{2})
        h_hsp = plot3(hsp(:,1),hsp(:,2),hsp(:,3),'.b');
        hold on
        h_fid2 = plot3(fid2(:,1),fid2(:,2),fid2(:,3),'sc');
        set(h_fid1,'MarkerFaceColor','c','MarkerSize',12,'MarkerEdgeColor','k');
    end
 
% cameramenu
zoom(1.5)
axis off;
axis equal;
rotate3d
drawnow

end
