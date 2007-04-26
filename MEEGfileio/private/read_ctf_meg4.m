function [meg] = read_ctf_meg4(fname, hdr, begsample, endsample, chanindx)

% READ_CTF_MEG4 reads specified samples from a CTF continous datafile
% It neglects all trial boundaries as if the data was acquired in
% non-continous mode.
%
% Use as
%   [meg] = read_ctf_meg4(filename, hdr, begsample, endsample, chanindx)
% where
%   filename	name of the datafile, including the .meg4 extension
%   header      with all data information (from read_ctf_meg4)
%   begsample   index of the first sample to read
%   endsample   index of the last sample to read
%   chanindx	index of channels to read (optional, default is all)
%
% See also READ_CTF_MEG4

% "VSM MedTech Ltd. authorizes the release into public domain under the
% GPL licence of the Matlab source code files "read_ctf_res4.m" and
% "read_ctf_meg4.m" by the authors of said files from the F.C. Donders
% Centre, Nijmegen, The Netherlands."

% Author(s): Jim McKay November 1999
% Last revision: Jim McKay
% Copyright (c) 1999-2000 CTF Systems Inc. All Rights Reserved.
%
% modifications Copyright (C) 2002, Ole Jensen
% modifications Copyright (C) 2003, Robert Oostenveld
%
% $Log: read_ctf_meg4.m,v $
% Revision 1.14  2005/10/04 15:52:09  roboos
% fixed bug that occured when data is split over more than two files (thanks to flodlan)
%
% Revision 1.13  2005/02/18 13:16:58  roboos
% VSM MedTech Ltd. authorised the release of this code in the public domain
% updated the copyrights, updated the help
%
% Revision 1.12  2004/06/21 19:33:08  roberto
% made 2GB warning dependent on global fb flag
%
% Revision 1.11  2003/07/23 15:02:27  roberto
% added check on valid input for read_ctf_meg4, other changes unknown
%
% Revision 1.10  2003/05/22 09:09:41  roberto
% fixed another bug for >2GB files when selected data in within one trial
%
% Revision 1.7  2003/05/21 13:52:29  roberto
% re-implemented support for >2GB files
% improved checking of input arguments
% fixed bug in chanindx indexing for raw data
%
% Revision 1.6  2003/05/19 15:18:50  roberto
% fixed bugs in memory-efficient reading of continuous data
%
% Revision 1.4  2003/04/17 12:37:41  roberto
% changed error for non-continuous files into warning
%
% Revision 1.3  2003/04/01 06:53:35  roberto
% added support for channel selection
% fixed bug with data allocation over multiple trials
%
% Revision 1.2  2003/03/27 08:30:54  roberto
% fixed bug in reading non-multiplexed trial data
% added error checking
%
% Revision 1.1  2003/03/26 13:34:05  roberto
% new implementation
%

% use global flag for feedback
global fb
if isempty(fb)
  fb = 0;
end

if begsample<1
  error('cannot read before the start of the data');
end

if endsample>hdr.nSamples*hdr.nChans*hdr.nTrials
  error('cannot read beyond the end of the data');
end

if begsample>endsample
  error('cannot read a negative number of samples');
end

if nargin<5
  % select all channels
  chanindx = 1:hdr.nChans;
end

if isempty(chanindx)
  error('no channels were specified for reading CTF data')
end

fid = fopen(fname,'r','ieee-be');

if fid == -1
  error('could not open datafile');
end

CTFformat = setstr(fread(fid,8,'char'))';
if (strcmp(CTFformat(1,1:7),'MEG41CP')==0),
    error('datafile is not in CTF MEG4 format')
end 

% the data is not channel multiplexed, but stored in trials
% in each trial first all samples for channel 1 are given, then all samples of channel 2 ...
% this means that we have to read each channel for each trial

begtrial = ceil(begsample/hdr.nSamples);
endtrial = ceil(endsample/hdr.nSamples);

% this counts the files and offset if the 2GB file size boundary is encountered
multifilecount  = 0;
multifileoffset = 0;
fseek(fid, 0, 'eof');
multifilelength = ftell(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read only the selected data, channel-wise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if endtrial==begtrial
  rawbegsample = begsample - (begtrial-1)*hdr.nSamples;
  rawendsample = endsample - (begtrial-1)*hdr.nSamples;
  for chan=1:length(chanindx)
    % jump to the begin of this channel in this trial
    channeloffset = 8*(multifilecount+1)+(begtrial-1)*hdr.nSamples*4*hdr.nChans+(chanindx(chan)-1)*hdr.nSamples*4-multifileoffset;
    if channeloffset>=multifilelength
      % data goes beyond 2GB file boundary, jump to the next file
      multifileoffset = multifileoffset + floor(channeloffset/multifilelength)*multifilelength;	% remember for the next trial
      multifilecount  = multifilecount +  floor(channeloffset/multifilelength);% increase the file counter
      channeloffset   = channeloffset  -  floor(channeloffset/multifilelength)*multifilelength+8;	% change the current offset, keep the header 
      nextname = sprintf('%s.%d_meg4', fname(1:(end-5)), multifilecount);
      if fb
        fprintf('data goes beyond 2GB file boundary, continuing with %s\n', nextname);
      end
      fclose(fid);
      fid = fopen(nextname,'r','ieee-be');
      fseek(fid, 0, 'eof');
      multifilelength = ftell(fid);				% determine the length of the current file
    end
    fseek(fid, channeloffset, 'bof');
    % jump to the first sample of interest
    fseek(fid, 4*(rawbegsample-1), 'cof');
    % read the data from this channel
    [tmp, count] = fread(fid,[(endsample-begsample+1),1],'int32');
    raw(:,chan) = tmp;
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read and concatenate the raw data of all trials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
  raw = zeros((endtrial-begtrial+1)*hdr.nSamples, length(chanindx));
  for trial=begtrial:endtrial

    if length(chanindx)==hdr.nChans & all(chanindx(:)'==1:hdr.nChans)
      % read the data from all channels
      rawbegsample = (trial-begtrial)*hdr.nSamples + 1;
      rawendsample = (trial-begtrial)*hdr.nSamples + hdr.nSamples;
      % jump to the begin of this trial
      trialoffset = 8*(multifilecount+1)+(trial-1)*hdr.nSamples*4*hdr.nChans-multifileoffset;
      if trialoffset>=multifilelength
        % data goes beyond 2GB file boundary, jump to the next file
        trialoffset   = trialoffset - multifilelength+8;	% change the current offset, keep the header 
        multifileoffset = multifileoffset + multifilelength;	% remember for the next trial
        multifilecount  = multifilecount + 1;			% increase the file counter
        nextname = sprintf('%s.%d_meg4', fname(1:(end-5)), multifilecount);
        if fb
          fprintf('data goes beyond 2GB file boundary, continuing with %s\n', nextname);
        end
        fclose(fid);
        fid = fopen(nextname,'r','ieee-be');
        fseek(fid, 0, 'eof');
        multifilelength = ftell(fid);				% determine the length of the current file
      end
      fseek(fid, trialoffset, 'bof');
      [tmp, count] = fread(fid,[hdr.nSamples,hdr.nChans],'int32');
      raw(rawbegsample:rawendsample,:) = tmp;
    else
      % read the data from the selected channels
      rawbegsample = (trial-begtrial)*hdr.nSamples + 1;
      rawendsample = (trial-begtrial)*hdr.nSamples + hdr.nSamples;
      for chan=1:length(chanindx)
        % jump to the begin of this channel in this trial
        channeloffset = 8*(multifilecount+1)+(trial-1)*hdr.nSamples*4*hdr.nChans+(chanindx(chan)-1)*hdr.nSamples*4-multifileoffset;
        if channeloffset>=multifilelength
          % data goes beyond 2GB file boundary, jump to the next file
          multifileoffset = multifileoffset + floor(channeloffset/multifilelength)*multifilelength; % remember for the next trial
          multifilecount  = multifilecount  + floor(channeloffset/multifilelength); % increase the file counter
          channeloffset   = channeloffset   - floor(channeloffset/multifilelength)*multifilelength+8;	% change the current offset, keep the header, multiply with the floor thing in the case of very big datasets 
          nextname = sprintf('%s.%d_meg4', fname(1:(end-5)), multifilecount);
          if fb
            fprintf('data goes beyond 2GB file boundary, continuing with %s\n', nextname);
          end
          fclose(fid);
          fid = fopen(nextname,'r','ieee-be');
          fseek(fid, 0, 'eof');
          multifilelength = ftell(fid);				% determine the length of the current file
        end
        fseek(fid, channeloffset, 'bof');
        [tmp, count] = fread(fid,[hdr.nSamples,1],'int32');
        raw(rawbegsample:rawendsample,chan) = tmp;
      end %chan

    end %read all channels
  end %trial

  % select the raw data corresponding to the samples of interest
  rawoffset    = (begtrial-1)*hdr.nSamples;
  rawbegsample = begsample - rawoffset;
  rawendsample = endsample - rawoffset;
  raw = raw(rawbegsample:rawendsample, :);
end
fclose(fid);

% multiply the dimensionless values with the calibration value
gain = hdr.gainV(chanindx);	% only for selected channels
meg = raw';			% transpose the raw data
for i=1:size(meg,1)
  meg(i,:) = gain(i)*meg(i,:);
end


