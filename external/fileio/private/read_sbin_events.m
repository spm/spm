function [EventCodes, segHdr, eventData] = read_sbin_events(filename)

% READ_SBIN_EVENTS reads the events information from an EGI segmented simple binary format file
%
% Use as
%   [EventCodes, segHdr, eventData] = read_sbin_events(filename)
% with
%   EventCodes      - if NEvent (from header_array) != 0, then array of 4-char event names
%   segHdr          - condition codes and time stamps for each segment
%   eventData       - if NEvent != 0 then event state for each sample, else 'none'
% and
%   filename    - the name of the data file
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
if bitand(version,1) == 0
	error('ERROR:  This is an unsegmented file.\n');
end;

if version > 7
	error('ERROR:  This is not a simple binary file.\n');
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

    if readNumSegments == 0
        readNumSegments = NSegments;            % If first and last segments not specified, read all of them.
    end;
    
	header_array 	= [version year month day hour minute second millisecond Samp_Rate NChan Gain Bits Range NumCategors, NSegments, NSamples, NEvent];

	for j = 1:firstSegment-1
        switch precision
            case 2
		        throwAway	= fread(NChan+NEvent, NSamples, 'int16');
            case 4
		        throwAway	= fread(NChan+NEvent, NSamples, 'single');
            case 6
		        throwAway	= fread(NChan+NEvent, NSamples, 'double');
        end
	end

	eventData	= zeros(NEvent,readNumSegments*NSamples);
    segHdr      = zeros(readNumSegments,2);
	
    if (NEvent ~= 0)
        for j = 1:readNumSegments
            [segHdr(j,1), count]	= fread(fid, 1,'int16');    %cell
            [segHdr(j,2), count]	= fread(fid, 1,'int32');    %time stamp
            switch precision
                case 2
                    [temp,count]	= fread(fid,[NChan+NEvent, NSamples],'int16');
                case 4
                    [temp,count]	= fread(fid,[NChan+NEvent, NSamples],'single');
                case 6
                    [temp,count]	= fread(fid,[NChan+NEvent, NSamples],'double');
            end
            eventData(:,((j-1)*NSamples+1):j*NSamples)	= temp( (NChan+1):(NChan+NEvent), 1:NSamples);
        end
    end
fclose(fid);
