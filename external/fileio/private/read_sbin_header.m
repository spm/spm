function [header_array, CateNames, CatLengths, preBaseline] = read_sbin_header(filename)

% READ_SBIN_HEADER reads the header information from an EGI segmented simple binary format file
%
% Use as
%   [header_array, CateNames, CatLengths, preBaseline] = read_sbin_header(filename)
% with
%   header_array     - differs between versions, read code for details
%   CateNames        - category names
%   CatLengths       - length of category names
%   preBaseline      - number of samples in the baseline prior to the baseline event
% and
%   filename    - the name of the data file
%
% Since there is no unique event code for the segmentation event, and hence the baseline period,
% the first event code in the list will be assumed to be the segmentation event.
% NetStation itself simply ignores possible baseline information when importing simple binary files.
%_______________________________________________________________________
%
%
% Modified from EGI's readEGLY.m with permission 2008-03-31 Joseph Dien
%

fid=fopen([filename],'r');
if fid==-1
    error('wrong filename')
end

version		= fread(fid,1,'int32');

if version > 7
    error('ERROR:  This is not a simple binary file.  Note that NetStation does not successfully directly convert EGIS files to simple binary format.\n');
end;

if bitand(version,1) == 0
    error('ERROR:  This is an unsegmented file.\n');
end;

precision = bitand(version,6);
if precision == 0
    error('File precision is not defined.');
end;

%		read header...
year		= fread(fid,1,'int16');
month		= fread(fid,1,'int16');
day			= fread(fid,1,'int16');
hour		= fread(fid,1,'int16');
minute		= fread(fid,1,'int16');
second		= fread(fid,1,'int16');
millisecond = fread(fid,1,'int32');
Samp_Rate	= fread(fid,1,'int16');
NChan		= fread(fid,1,'int16');
Gain 		= fread(fid,1,'int16');
Bits 		= fread(fid,1,'int16');
Range 		= fread(fid,1,'int16');
NumCategors	= fread(fid,1,'int16');
for j = 1:NumCategors
    CatLengths(j)	= fread(fid,1,'int8');
    for i = 1:CatLengths(j)
        CateNames(j,i)	= char(fread(fid,1,'char'));
    end
end
NSegments	= fread(fid,1,'int16');
NSamples	= fread(fid,1,'int32');			% samples per segment
NEvent		= fread(fid,1,'int16');			% num events per segment
EventCodes = [];
for j = 1:NEvent
    EventCodes(j,1:4)	= char(fread(fid,[1,4],'char'));
end

header_array 	= [version year month day hour minute second millisecond Samp_Rate NChan Gain Bits Range NumCategors, NSegments, NSamples, NEvent];

preBaseline=0;
if NEvent > 0

    %read the first segment to determine baseline length (assuming all segments have the same baseline).
    j=1;
    [temp1]	= fread(fid, 1,'int16');    %cell
    [temp2]	= fread(fid, 1,'int32');    %time stamp
    switch precision
        case 2
            [temp,count]	= fread(fid,[NChan+NEvent, NSamples],'int16');
        case 4
            [temp,count]	= fread(fid,[NChan+NEvent, NSamples],'single');
        case 6
            [temp,count]	= fread(fid,[NChan+NEvent, NSamples],'double');
    end
    eventData = temp( (NChan+1):(NChan+NEvent), 1:NSamples);
    theEvent=find(eventData(1,:)>0); %assume the first event code is the segmentation event
    theEvent=theEvent(1);
    preBaseline=theEvent-1;
    if preBaseline == -1
        preBaseline =0;
    end; 
end

fclose(fid);