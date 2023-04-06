function [D,L] = spm_opm_create(S)
% Read magnetometer data and optionally set up forward model
% FORMAT [D,L] = spm_opm_create(S)
%   S               - input structure
% Optional fields of S:
% SENSOR LEVEL INFO
%   S.data          - filepath/matrix(nchannels x timepoints)  - Default:required
%   S.channels      - channels.tsv file                        - Default: REQUIRED
%   S.fs            - Sampling frequency (Hz)                  - Default: REQUIRED if S.meg is empty
%   S.meg           - meg.json file                            - Default: REQUIRED if S.fs is empty
%   S.precision     - 'single' or 'double'                     - Default: 'single'
% SOURCE LEVEL INFO
%   S.coordsystem   - coordsystem.json file                    - Default: transform between sensor space and anatomy is identity
%   S.positions     - positions.tsv file                       - Default: no Default
%   S.sMRI          - Filepath to  MRI file                    - Default: no Default
%   S.template      - Use SPM canonical template               - Default: 0
%   S.headhape      - .pos file for better template fit        - Default:
%   S.cortex        - Custom cortical mesh                     - Default: Use inverse normalised cortical mesh
%   S.scalp         - Custom scalp mesh                        - Default: Use inverse normalised scalp mesh
%   S.oskull        - Custom outer skull mesh                  - Default: Use inverse normalised outer skull mesh
%   S.iskull        - Custom inner skull mesh                  - Default: Use inverse normalised inner skull mesh
%   S.voltype       - Volume conducter Model type              - Default: 'Single Shell'
%   S.meshres       - mesh resolution(1,2,3)                   - Default: 1
%   S.lead          - flag to compute lead field               - Default: 0
% Output:
%   D               - MEEG object (also written to disk)
%   L               - Lead field (also written on disk)
%__________________________________________________________________________

% Tim Tierney
% Copyright (C) 2018-2023 Wellcome Centre for Human Neuroimaging


spm('FnBanner', mfilename);

%-Set default values
%--------------------------------------------------------------------------
if ~isfield(S, 'voltype'),     S.voltype = 'Single Shell';  end
if ~isfield(S, 'meshres'),     S.meshres = 1;               end
if ~isfield(S, 'sMRI'),        S.sMRI = [];                 end
if ~isfield(S, 'scalp'),       S.scalp = [];                end
if ~isfield(S, 'template'),    S.template = 0;              end
if ~isfield(S, 'cortex'),      S.cortex = [];               end
if ~isfield(S, 'iskull'),      S.iskull = [];               end
if ~isfield(S, 'oskull'),      S.oskull = [];               end
if ~isfield(S, 'fname'),       S.fname = 'sim_opm';         end
if ~isfield(S, 'path'),        S.path = [];                 end
if ~isfield(S, 'precision'),   S.precision = 'single';      end
if ~isfield(S, 'lead'),        S.lead = 0;                  end
if ~isfield(S, 'headshape'),   S.headshape = [];            end
if ~isfield(S, 'coordsystem'), S.coordsystem = [];          end

%- identify Binary File
%----------------------------------------------------------------------
try % work out if data is a matrix or a file
    [direc, dataFile] = fileparts(S.data);
    binData = true;
    
    %- Cerca magnetics support
    %----------------------------------------------------------------------
    if strcmpi(num2str(S.data(end-4:end)),'.cMEG')
        forward = 0;
        subjectSource=0;
        subjectNoStruct = 0;
        template = 0;
        warning(['Cerca magnetics data currently requires separate'...
            ' coregistration and forward modelling.']);
        args = [];
        args.filename = S.data;
        if isfield(S, 'positions')
            positions = 1;
            args.positions = S.positions;
        else
            positions = 0;
        end
        Data = read_cMEG_data(args);
        S.channels = Data.channels;
        if isfield(Data,'position')
            S.positions = Data.position;
        end
        S.fs = Data.meg.SamplingFrequency;
        binData = false;
        S.data = Data.data;
    end
    
catch % if not readable check if it is numeric
    if ~isa(S.data,'numeric')
        error('A valid dataset or file was not supplied.');
    end
    binData = false;
    direc = pwd();
    dataFile = S.fname;
    if ~isfield(S, 'channels')
        error('A channels.tsv file must be supplied.');
    end
end

%- identify potential BIDS Files
%--------------------------------------------------------------------------
base = strsplit(dataFile,'meg');
chanFile = fullfile(direc,[base{1},'channels.tsv']);
megFile = fullfile(direc,[base{1},'meg.json']);
posFile = spm_select('FPList',direc,[base{1},'positions.tsv']);
if isempty(S.coordsystem)
    S.coordsystem = fullfile(direc,[base{1},'coordsystem.json']);
end

%- Check for channel Info
%--------------------------------------------------------------------------

chanfieldexists = isfield(S,'channels');

if chanfieldexists   
    ischanstruct = isstruct(S.channels);
end

if (chanfieldexists && ischanstruct)
    % use channel struct if it exists
    channels = S.channels;
else
    try % to load a channels file
        channels = spm_load(S.channels);
    catch
        try  % use channel struct if supplied
            channels = spm_load(chanFile);
        catch % create channel struct
            error('A valid channels.tsv file or struct was not found.');
        end
    end
end

%- Check for MEG Info
%--------------------------------------------------------------------------
try % to load a meg file
    meg = spm_load(S.meg);
catch
    try % to load a BIDS meg file
        meg = spm_load(megFile);
    catch
        try % to use meg struct
            meg = S.meg;
        catch
            try % to use S.fs argument to get sampling frequency
                meg = [];
                meg.SamplingFrequency = S.fs;
            catch
                error('A meg.json file is required if S.fs is empty.');
            end
        end
    end
end

%- Position File check
%--------------------------------------------------------------------------
try % to load a channels file
    posOri = spm_load(S.positions);
    positions = 1;
catch
    try % to load a BIDS channel file
        posOri = spm_load(posFile);
        positions = 1;
    catch
        try % to assign a BIDS struct of positions
            if (isstruct(S.positions))
                posOri = S.positions;
                positions = 1;
            else
                positions = 0;
            end
        catch
            warning('No sensor position information found.');
            positions = 0;
        end
    end
end

%- work out data size
%--------------------------------------------------------------------------
nChans = size(channels.name,1);

if binData
    fprops = dir(S.data);
    if strcmp(S.precision,'single')
        bytesPerSample = 4;
    else
        bytesPerSample = 8;
    end
    nSamples = fprops.bytes/(nChans*bytesPerSample);
    nTrials = 1;
else
    nSamples = size(S.data,2);
    nTrials = size(S.data,3);
end

%- Check for a custom save path
%--------------------------------------------------------------------------
if ~isempty(S.path)
    if exist(S.path,'dir')
        direc = S.path;
    else
        error('Specified output directory does not exist.');
    end
end % will use original direc variable otherwise

%- Create SPM object
%--------------------------------------------------------------------------
D = meeg(nChans,nSamples,nTrials);
D = fsample(D,meg.SamplingFrequency);
D = fname(D,[dataFile,'.mat']);
D = path(D,direc);
D = chanlabels(D,1:size(D,1),channels.name);
D = units(D,1:size(D,1),channels.units);
D = chantype(D,1:size(D,1),channels.type);

%- Overwrite and Save
%--------------------------------------------------------------------------
ma = fullfile(direc,[dataFile,'.mat']);
da = fullfile(direc,[dataFile,'.dat']);

if exist(fullfile(path(D),fname(D)),'file') == 2
    delete(ma);
    delete(da);
end
D.save();

%- Reformat data according to channel info
%--------------------------------------------------------------------------
D = blank(D,[dataFile,'.dat']);

if binData
    maxMem = 100e6;
    samplesPerChunk = round(maxMem/(8*nChans));
    begs = 1:samplesPerChunk:nSamples;
    ends = begs+samplesPerChunk-1;
    
    if ends(end) > nSamples
        ends(end) = nSamples;
    end
    
    dat = fopen(S.data);
    for j=1:length(begs)
        asamplesPerChunk = length(begs(j):ends(j));
        inds = begs(j):ends(j);
        D(:,inds,1) = fread(dat,[nChans,asamplesPerChunk],S.precision,0,'b');
    end
    fclose(dat);
else
    % insert provided data
    D(1:nChans,1:nSamples,1:nTrials) = S.data;
    D.save();
end

%-Place sensors in object
%--------------------------------------------------------------------------
if positions
    pos = [posOri.Px,posOri.Py,posOri.Pz];
    ori = [posOri.Ox,posOri.Oy,posOri.Oz];
    cl = posOri.name;
    
    [sel1, sel2] = match_str(cl,channels.name);
    
    grad = [];
    grad.label = cl;
    grad.coilpos = pos;
    grad.coilori = ori;
    grad.tra = eye(numel(grad.label));
    grad.chanunit = repmat({'T'}, numel(grad.label), 1);
    grad.chantype = lower({channels.type{sel2}}');
    grad = ft_datatype_sens(grad, 'amplitude', 'T', 'distance', 'mm');
    D = sensors(D, 'MEG', grad);
    save(D);
    
    
    %- 2D view based on mean orientation of sensors
    n1 = mean(grad.coilori); n1= n1./sqrt(dot(n1,n1));
    t1 = cross(n1,[0 0 1]);
    t2 = cross(t1,n1);
    pos2d = zeros(size(grad.coilpos,1),2);
    for i=1:size(grad.coilpos,1)
        pos2d(i,1) = dot(grad.coilpos(i,:),t1);
        pos2d(i,2) = dot(grad.coilpos(i,:),t2);
    end
    
    nMEG = length(indchantype(D,'MEG'));
    if nMEG ~= size(pos2d,1)
        m1 = '2D positions could not be set as there are ';
        m2 = num2str(nMEG);
        m3 = ' channels but only ';
        m4 = num2str(size(pos2d,1));
        m5 = ' channels with position information.';
        message = [m1,m2,m3,m4,m5];
        warning(message);
    else
        args = [];
        try
            args.xy = Data.xy';
            args.label = Data.xylabs;
        catch
            args.xy = pos2d';
            args.label = grad.label;
        end
        args.D = D;
        args.task = 'setcoor2d';
        D = spm_eeg_prep(args);
        D.save;
    end
end

%- prep for headmodel processing
%--------------------------------------------------------------------------
if positions
    
    % auto-populate fields for spm_opm_headmodel
    targets = {'template','sMRI','coordsystem','headshape',...
        'meshres','iskull','oskull','scalp','cortex','voltype','lead'};
    args = struct;
    args.D = D;
    for ii = 1:numel(targets)
        args.(targets{ii}) = S.(targets{ii});
    end
    
    if S.lead
        [D,L] = spm_opm_headmodel(args);
    else
        D = spm_opm_headmodel(args);
    end
    D.save;
end
D.save();

fprintf('%-40s: %30s\n','Completed',spm('time'));

end


function [bids] = read_cMEG_data(S)
% Read OPM data from cerca magnetics in bids structs suitable for
% conversion to SPM format
% FORMAT bids = read_cMEG_data(S)
%   S               - input structure
% Optional fields of S:
% SENSOR LEVEL INFO
%   S.filename      - filepath                          - Default:required
%   S.positions     - cerca positions.mat file          - Default: not used
% Output:
%  bids             - struct containing info for channels.tsv,meg.json and
%                     positions.tsv
%__________________________________________________________________________


if strcmpi(num2str(S.filename(end-4:end)),'.cMEG')
    filename = S.filename(1:end-5);
end

if ~isfield(S, 'positions')
    pos = false;
else
    pos = true;
end

%disp('Get data')
fname = [filename '.cMEG'];
fid = fopen(fname,'rb','ieee-be');
finfo = dir(fname);
fsize = finfo.bytes;
Adim_conv = [2^32; 2^16; 2^8; 1]; % Array dimension conversion table

%disp('Preallocation')
I = 0;
while fsize ~= ftell(fid)
    dims = [];
    for n = 1:2 % 2 dimensions (Nch x Time)
        temp = fread(fid,4);
        temp = sum(Adim_conv.*temp);
        dims = [dims,temp];
    end
    I = I + 1;
    fread(fid,prod(dims),'double',0,'ieee-be');  % Skip the actual data (much faster for some reason)
end
fseek(fid,0,'bof'); % Reset the cursor
data1 = repmat({NaN*ones(dims)},I,1);  % Preallocate space, assumes each section is the same

%disp('Load and parse data')
for j = 1:I
    dims = [];  % Actual array dimensions
    for n = 1:2 % 2 dimensions (Nch x Time)
        temp = fread(fid,4);
        temp = sum(Adim_conv.*temp);
        dims = [dims,temp];
    end
    clear temp n
    temp = fread(fid,prod(dims),'double',0,'ieee-be');
    
    % Reshape the data into the correct array configuration
    for n = 1:2
        temp = reshape(temp,dims(2-n+1),[])';
    end
    data1{j} = temp;  % Save the data
end
fclose(fid);  % Close the file

%disp('Reshape into sensible order')
data = NaN*zeros(size(data1{1},2),size(data1{1},1).*size(data1,1));
for n = 1:size(data1,1)
    %clc
    %     disp(100*n/length(data1))
    data_range = [(n-1)*size(data1{1},1)+1:n*size(data1{1},1)];
    data(:,data_range) = [data1{n}']; % data in Nch+triggers+1 x time, but channel one is the time
end

%disp('Read session info')
Session_info_file = [filename(1:15) '_SessionInfo.txt'];
fidSI = fopen(Session_info_file);
finfoSI = dir(Session_info_file);
fsizeSI = finfoSI.bytes;
count = 0;
while fsizeSI ~= ftell(fidSI)
    count = count + 1;
    Session_info{count,1} = fgetl(fidSI);
end
fclose(fidSI);

channels = spm_load([filename(1:15) '_channels.tsv']);
Chan_names = channels.name;
for i = 1:size(Chan_names,1)
basespl = strsplit(channels.name{i},'[');
Chan_names{i}= [basespl{1} basespl{2}(1)];
end
channels.name=Chan_names;

f = round(1/(data(1,2)-data(1,1)));
time = linspace(0,size(data,2)./f,size(data,2));
data = data(2:end,:);

% Convert to fT (1x Mode, so 2.7 V/nT)
gainIndex = find(not(cellfun('isempty',strfind(Session_info,'OPM Gain'))));
gain = Session_info{gainIndex};
gainString = strsplit(gain, ': ');
gainFac = str2num((gainString{2}(1:4)));
opminds = strcmp(channels.type,'MEGMAG');
data(opminds,:) = (1e6.*data(opminds,:))./(2.7*gainFac);
channels.units(opminds)={'fT'};

%- positions
%--------------------------------------------------------------------------
if pos
    positions=[];
    helm = spm_load(S.positions);
    
    [~,b] = spm_match_str(Chan_names,helm.Sensor);
    positions.name= helm.Sensor(b);
    
    if(iscell(helm.Px(b)))
        positions.Px = str2double(helm.Px(b));
    else
        positions.Px =helm.Px(b);
    end
    
    if(iscell(helm.Py(b)))
        positions.Py = str2double(helm.Py(b));
    else
        positions.Py =helm.Py(b);
    end
    
    if(iscell(helm.Pz(b)))
        positions.Pz = str2double(helm.Pz(b));
    else
        positions.Pz =helm.Pz(b);
    end
    
    if(iscell(helm.Ox(b)))
        positions.Ox = str2double(helm.Ox(b));
    else
        positions.Ox =helm.Ox(b);
    end
     if(iscell(helm.Oy(b)))
        positions.Oy = str2double(helm.Oy(b));
    else
        positions.Oy =helm.Oy(b);
     end
     
     if(iscell(helm.Oz(b)))
        positions.Oz = str2double(helm.Oz(b));
    else
        positions.Oz =helm.Oz(b);
    end
    
    xy = zeros(size(positions.name,1),2);
    xylabs = positions.name;
    xy(:,1) = helm.Layx(b);
    xy(:,2) = helm.Layy(b);
    
    bids.xy = xy;
    bids.xylabs = xylabs;
    bids.position = positions;
end


%- MEG.json - only crucial fields populated
%--------------------------------------------------------------------------
meg = [];
meg.SamplingFrequency = f;

%- 1 output struct
%--------------------------------------------------------------------------
bids.data = data;
bids.channels = channels;
bids.meg = meg;

end
