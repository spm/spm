function T = spm_img_type(P)
% figure out image type
% FORMAT T = spm_img_type(P)
% P - filename
% T - format code
%____________________________________________________________________________
%
% For filename P the file format is divined and returned.  
% Format codes have the following meanings:
%
%  0 = Unknown
%  1 = Analyze
%  2 = Bad Analyze
%  3 = CTI
%  4 = Bad CTI  
%  7 = AIR Parameter file
%  8 = Bad AIR file
%
% from spm_img_type.m, v1.2, UPMC PET Modified SPM (4/26/95)
% %W% %E%

% These correspond to return values of which_imgtype() in pAIR
UNKNOWN		= 0;
ANALYZE		= 1;
BAD_ANALYZE	= 2;
CTI		= 3;
BAD_CTI		= 4;
AIR_PARAM	= 7;
BAD_AIR_PARAM	= 8;


Pb = spm_clip_suffix(P);

%
% See if it's Analyze
%
% If correct length header exists (and header length field is set
% correctly) AND and image file exists, a file is type ANALYZE.
%
% If header file exists and its the wrong length OR if the image file
% doesn't exist OR if the sizeof_hdr field is set wrong, the file 
% type is BAD_ANALYZE
%
AnalHdrSz = 348;

PA    = [Pb '.hdr'];
fid   = fopen(PA);
if (fid ~= -1)
	T = BAD_ANALYZE;
	% Can we seek?
	if fseek(fid,0,'eof') ~= -1
		% Is the header the right length?
		if ftell(fid) == AnalHdrSz
			% Is sizeof_hdr field set right?
			fseek(fid,0,'bof');
			if fread(fid,1,'long') == AnalHdrSz
				% Does the .img file exist?
				fid2 = fopen([Pb '.img']);
				if fid2 ~= -1
					T = ANALYZE;
					fclose(fid2);
				end
			end
		end
		
	end
	fclose(fid);
	return
end


%
% See if it's an CTI file
%
% A CTI file must be a multiple of CTIblkSz AND have the specified
% Magic Numbers at the specified offsets.  The most reliable magic numbers
% I could come up with were software version and datatype in the main
% header, but alas, these will be changing with ECAT 7.0, so you can
% specify multple valid values for a given offset.
%
% I can think of no easy way to detect a BAD_CTI file (if files lengths
% of multiple 512 were less common, I'd use that)
%
CTIblkSz = 512;
% at offset 48 is software version (could be 6 or 7)
% at offset 50 is datatype (could be 2 or 6 {VAX or SUN I2})
% at offset 54 is filetype (must be 2 {IMAGE file})
CTImagicN = [	50 2 6;  
		54 2 2 ];  

if (strcmp(P, Pb))
	% If there was no suffix, try '.img'
	PA = [Pb '.img'];
else
	PA = P;
end
fid   = fopen(PA);
if (fid ~= -1)

	% Can we seek?
	if fseek(fid,0,'eof') ~= -1

	    % Is the file a unit of CTIblkSz?
	    if rem(ftell(fid),CTIblkSz) == 0

		% Get values at all offsets
		MgN = [];
		nMgN = size(CTImagicN,1);
	    	for i  = 1:nMgN
	    	    if fseek(fid,CTImagicN(i,1),'bof') ~= -1
	                byt = fread(fid,2,'uchar');
			MgN = [MgN; byt(1)+byt(2)*255];
	            end	                
	    	end

	        % Don't need the offsets anymore
		d = size(CTImagicN,2);
                CTImagicN = CTImagicN(:,2:d);

		if all(any(MgN(:,ones(1,d-1))' == CTImagicN'))
				T = CTI;
		end		

	    end

	end
	fclose(fid);
	if (T == CTI) | (T == BAD_CTI)
		return
	end
end


%
% See if it's an AIR file
%
%
AirParamSzs = [	688;	% Second AIR file size is from a old UPMC beta AIR.
              	808];	% Word from Woods: The next AIR will have yet 
			% a different file size

PA    = [Pb '.air'];
fid   = fopen(PA);
if (fid ~= -1)
	T = BAD_AIR_PARAM;
 	if fseek(fid,0,'eof') ~= -1
		Sz = ftell(fid);
		if Sz > -1;
			d = size(AirParamSzs,1);
			if any(Sz(ones(1,d),:) == AirParamSzs)
				T = AIR_PARAM;
			end
		end
	end
	fclose(fid);
	return
end


T = UNKNOWN;

return
