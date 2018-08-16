function M = mz3_read(filename)
% Read MZ3-formatted data from disk
% FORMAT M = mz3_read(filename)
%
% filename - MZ3-formatted file name
% M        - data structure
%__________________________________________________________________________
% 
% MZ3 Format Specification:
% https://github.com/neurolabusc/surf-ice/tree/master/mz3
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: mz3_read.m 7399 2018-08-16 12:04:14Z guillaume $


M = struct();

%-Open file
%--------------------------------------------------------------------------
fid = fopen(filename,'r','ieee-le');
if fid == -1
    error('Cannot open %s.',filename);
end

%-Load (uncompressed) data
%--------------------------------------------------------------------------
magic = fread(fid,[1 2],'uint8');
if isequal(magic,[77 90]) % uncompressed 0x4D5A, 23117, 'MZ'
    frewind(fid);
    D = fread(fid,Inf,'*uint8');
elseif isequal(magic,[31 139]) % GZipped compressed 0x1F8B
    % Read gzip file
    % http://zlib.org/rfc-gzip.html
    member = fread(fid,8,'uint8'); % |CM |FLG|     MTIME     |XFL|OS |
    if member(1) ~= 8
        fclose(fid);
        error('Unknown compression method.');
    end
    flag = uint8(member(2));
    if bitget(flag,2) % (if FLG.FEXTRA set)
        xlen = fread(fid,1,'uint16'); % | XLEN  |
        fread(fid,xlen,'uint8'); % |...XLEN bytes of "extra field"...| 
    end
    if bitget(flag,3) % (if FLG.FNAME set)
        while true
            if fread(fid,1,'uint8') == 0
                break; % |...original file name, zero-terminated...|
            end
        end
    end
    if bitget(flag,4) % (if FLG.FCOMMENT set)
        while true
            if fread(fid,1,'uint8') == 0
                break; % |...file comment, zero-terminated...|
            end
        end
    end
    if bitget(flag,1) % (if FLG.FHCRC set)
        fread(fid,2,'uint8'); % | CRC16 |
    end
    % |...compressed blocks...|     CRC32     |     ISIZE     |
    d = fread(fid,Inf,'*uint8');
    isize = typecast(d(end-3:end),'uint32');
    D = zstream('d',d);
    if numel(D) ~= isize % size of the uncompressed input data modulo 2^32
        warning('Decompression failed.');
    end
else
    fclose(fid);
    error('Invalid file %s.',filename);
end

%-Close file
%--------------------------------------------------------------------------
fclose(fid);


%-Parse data
%==========================================================================
attr  = typecast(D(3:4),'uint16');
Nf    = typecast(D(5:8),'uint32');
Nv    = typecast(D(9:12),'uint32');
Nskip = typecast(D(13:16),'uint32') + 16; % total header size

if bitget(attr,1) % isFace
    M.faces = reshape(typecast(D(Nskip+(1:Nf*3*4)),'uint32'),3,[])' + 1;
    Nskip = Nskip + Nf*3*4;
end
if bitget(attr,2) % isVert
    M.vertices = reshape(typecast(D(Nskip+(1:Nv*3*4)),'single'),3,[])';
    Nskip = Nskip + Nv*3*4;
end
if bitget(attr,3) % isRGBA
    M.cdata = reshape(D(Nskip+(1:Nv*4)),4,[])';
    Nskip = Nskip + Nv*4;
end
if bitget(attr,4) % isScalar
    M.cdata = typecast(D(Nskip+(1:Nv*4)),'single');
    Nskip = Nskip + Nv*4;
end
