function [dat] = read_neuralynx_sdma(dataset, begsample, endsample, chanindx);

% READ_NEURALYNX_SDMA read specified channels and samples from a Neuralynx splitted DMA dataset
%
% Use as
%    [hdr] = read_neuralynx_sdma(dataset)
%    [dat] = read_neuralynx_sdma(dataset, begsample, endsample, chanindx)
%
% The splitted DMA dataset is not a formal Neuralynx format, but at
% the FCDC we use it in conjunction with SPIKEDOWNSAMPLE. The dataset
% directory contains files, one for each channel, each containing a
% 8-byte header followed by 32-bit integer values for all samples.

% Copyright (C) 2006, Robert Oostenveld
%
% $Log: read_neuralynx_sdma.m,v $
% Revision 1.5  2007/02/21 09:57:13  roboos
% fixed bug in TimeStampsPerSample (should be divided by nsamples-1)
% implemented re-ordering of channel names (i.e. not on ascii-order) to match the channel ordering to the original unspiltted DMA file
%
% Revision 1.4  2007/01/09 09:42:29  roboos
% read low and high timestamp as uint32, use subfunction to combine them into a uint64
%
% Revision 1.3  2006/12/12 11:39:57  roboos
% cleaned up code, changed handling of channel names (now based on file extension), directory listing determines channel ordering
%
% Revision 1.2  2006/12/04 20:26:02  roboos
% change in whitespace
%
% Revision 1.1  2006/04/24 12:21:14  roboos
% new implementation
%

needhdr = (nargin==1);
needdat = (nargin>=2);

if ~isdir(dataset)
  error('dataset should be a directory');
end

% determine the content of the dataset directory
l = dir(dataset);

% determine the correspondence between files and channels
label = {};
[pd, nd, xd] = fileparts(dataset);
for i=1:length(l)
  [pf, nf, xf] = fileparts(l(i).name);
  % compare the name of the channel with the name of the dataset
  if strcmp(nf, nd) && ~isempty(xf)
    label{i}   = xf(2:end);
    filesel(i) = 1;
  else
    label{i}   = [];
    filesel(i) = 0;
  end
end
label = label(find(filesel));
label = label(:);

% there is a preferred ordering, which corresponds with the channel order
% in the non-splitted DMA file

preferred = {
  'stx'
  'pid'
  'siz'
  'tsh'
  'tsl'
  'cpu'
  'ttl'
  'x01'
  'x02'
  'x03'
  'x04'
  'x05'
  'x06'
  'x07'
  'x08'
  'x09'
  'x10'
  'csc001'
  'csc002'
  'csc003'
  'csc004'
  'csc005'
  'csc006'
  'csc007'
  'csc008'
  'csc009'
  'csc010'
  'csc011'
  'csc012'
  'csc013'
  'csc014'
  'csc015'
  'csc016'
  'csc017'
  'csc018'
  'csc019'
  'csc020'
  'csc021'
  'csc022'
  'csc023'
  'csc024'
  'csc025'
  'csc026'
  'csc027'
  'csc028'
  'csc029'
  'csc030'
  'csc031'
  'csc032'
  'csc033'
  'csc034'
  'csc035'
  'csc036'
  'csc037'
  'csc038'
  'csc039'
  'csc040'
  'csc041'
  'csc042'
  'csc043'
  'csc044'
  'csc045'
  'csc046'
  'csc047'
  'csc048'
  'csc049'
  'csc050'
  'csc051'
  'csc052'
  'csc053'
  'csc054'
  'csc055'
  'csc056'
  'csc057'
  'csc058'
  'csc059'
  'csc060'
  'csc061'
  'csc062'
  'csc063'
  'csc064'
  'csc065'
  'csc066'
  'csc067'
  'csc068'
  'csc069'
  'csc070'
  'csc071'
  'csc072'
  'csc073'
  'csc074'
  'csc075'
  'csc076'
  'csc077'
  'csc078'
  'csc079'
  'csc080'
  'csc081'
  'csc082'
  'csc083'
  'csc084'
  'csc085'
  'csc086'
  'csc087'
  'csc088'
  'csc089'
  'csc090'
  'csc091'
  'csc092'
  'csc093'
  'csc094'
  'csc095'
  'csc096'
  'csc097'
  'csc098'
  'csc099'
  'csc100'
  'csc101'
  'csc102'
  'csc103'
  'csc104'
  'csc105'
  'csc106'
  'csc107'
  'csc108'
  'csc109'
  'csc110'
  'csc111'
  'csc112'
  'csc113'
  'csc114'
  'csc115'
  'csc116'
  'csc117'
  'csc118'
  'csc119'
  'csc120'
  'csc121'
  'csc122'
  'csc123'
  'csc124'
  'csc125'
  'csc126'
  'csc127'
  'csc128'
  'csc129'
  'csc130'
  'csc131'
  'csc132'
  'csc133'
  'csc134'
  'csc135'
  'csc136'
  'csc137'
  'csc138'
  'csc139'
  'csc140'
  'csc141'
  'csc142'
  'csc143'
  'csc144'
  'csc145'
  'csc146'
  'csc147'
  'csc148'
  'csc149'
  'csc150'
  'csc151'
  'csc152'
  'csc153'
  'csc154'
  'csc155'
  'csc156'
  'csc157'
  'csc158'
  'csc159'
  'csc160'
  'csc161'
  'csc162'
  'csc163'
  'csc164'
  'csc165'
  'csc166'
  'csc167'
  'csc168'
  'csc169'
  'csc170'
  'csc171'
  'csc172'
  'csc173'
  'csc174'
  'csc175'
  'csc176'
  'csc177'
  'csc178'
  'csc179'
  'csc180'
  'csc181'
  'csc182'
  'csc183'
  'csc184'
  'csc185'
  'csc186'
  'csc187'
  'csc188'
  'csc189'
  'csc190'
  'csc191'
  'csc192'
  'csc193'
  'csc194'
  'csc195'
  'csc196'
  'csc197'
  'csc198'
  'csc199'
  'csc200'
  'csc201'
  'csc202'
  'csc203'
  'csc204'
  'csc205'
  'csc206'
  'csc207'
  'csc208'
  'csc209'
  'csc210'
  'csc211'
  'csc212'
  'csc213'
  'csc214'
  'csc215'
  'csc216'
  'csc217'
  'csc218'
  'csc219'
  'csc220'
  'csc221'
  'csc222'
  'csc223'
  'csc224'
  'csc225'
  'csc226'
  'csc227'
  'csc228'
  'csc229'
  'csc230'
  'csc231'
  'csc232'
  'csc233'
  'csc234'
  'csc235'
  'csc236'
  'csc237'
  'csc238'
  'csc239'
  'csc240'
  'csc241'
  'csc242'
  'csc243'
  'csc244'
  'csc245'
  'csc246'
  'csc247'
  'csc248'
  'csc249'
  'csc250'
  'csc251'
  'csc252'
  'csc253'
  'csc254'
  'csc255'
  'csc256'
  'crc'
  };

[sel1, sel2] = match_str(preferred, label);
if length(sel2)==length(label)
  % all reorder the labels/files
  label = label(sel2);
else
  % not all files in this SDMA dataset could be accounted for in the
  % preference list, hence do not attempt to apply the preferred ordering
end

if needhdr
  % some parts of the header has to be hardcoded, since this dataset does not contain header information
  hdr.Fs           = 32556;     % sampling frequency
  hdr.nSamplesPre  = 0;         % number of pre-trigger samples in each trial
  hdr.nTrials      = 1;         % number of trials
  hdr.label        = label;
  hdr.nChans       = length(label);

  % determine the number of samples by looking at the file sizes in the directory
  s = median(cell2mat({l(find(filesel)).bytes}));
  hdr.nSamples = (s-8)/4;

  % determine the first and last timestamp, by reading them from the timestamp channels
  [p, f, x] = fileparts(dataset);
  fid = fopen(fullfile(dataset, [f, '.tsl']), 'rb', 'ieee-le');
  fseek(fid, 8, 'cof');
  beg_tsl = fread(fid, 1, 'uint32=>uint32');
  fseek(fid, 8+(hdr.nSamples-1)*4, 'bof');
  end_tsl = fread(fid, 1, 'uint32=>uint32');
  fclose(fid);
  fid = fopen(fullfile(dataset, [f, '.tsh']), 'rb', 'ieee-le');
  fseek(fid, 8, 'cof');
  beg_tsh = fread(fid, 1, 'uint32=>uint32');
  fseek(fid, 8+(hdr.nSamples-1)*4, 'bof');
  end_tsh = fread(fid, 1, 'uint32=>uint32');
  fclose(fid);
  hdr.FirstTimeStamp = timestamp_neuralynx(beg_tsl, beg_tsh);
  hdr.LastTimeStamp  = timestamp_neuralynx(end_tsl, end_tsh);
  hdr.TimeStampPerSample = double(hdr.LastTimeStamp-hdr.FirstTimeStamp)./(hdr.nSamples-1);  % this should be double, since it can be fractional

  % only return the header information
  dat = hdr;

else
  % allocate memory for all data
  dat = zeros(length(chanindx), endsample-begsample+1);

  % read all channels, one small chunk at at time, and write it to seperate files
  [p, f, x] = fileparts(dataset);
  for i=1:length(chanindx)
    % the name of the file can be reconstructed from the dataset and the channel label
    j = chanindx(i);
    chanfile = fullfile(dataset, [f '.' label{j}]);
    fid = fopen(chanfile, 'rb', 'ieee-le');
    % skip the 8-byte channel header
    fseek(fid, 8+begsample-1, 'cof');
    dat(i,:) = fread(fid, endsample-begsample+1, 'int32');
    fclose(fid);
  end
end
