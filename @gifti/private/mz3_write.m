function mz3_write(M,filename,fmt)
% Write MZ3-formatted data from disk
% FORMAT mz3_write(M,filename,fmt)
%
% M        - data structure
% filename - MZ3 output filename [Default: 'untitled']
% fmt      - compress data [Default: false]
%__________________________________________________________________________
% 
% MZ3 Format Specification:
% https://github.com/neurolabusc/surf-ice/tree/master/mz3
%__________________________________________________________________________
% Copyright (C) 2018 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: mz3_write.m 7470 2018-11-01 17:40:18Z guillaume $


if nargin < 3, fmt = false; end

%-Open file
%--------------------------------------------------------------------------
fid = fopen(filename,'wb','ieee-le');
if fid == -1
    error('Cannot open %s.',filename);
end

%-Write header
%--------------------------------------------------------------------------
attr = uint16(0);
Nf = 0; Nv = 0; Ns = 0;
if isfield(M,'faces')
    attr = bitset(attr,1);
    Nf = size(M.faces,1);
end
if isfield(M,'vertices')
    attr = bitset(attr,2);
    Nv = size(M.vertices,1);
end
if isfield(M,'cdata')
    if isa(M.cdata,'uint8') && size(M.cdata,2) == 4
        attr = bitset(attr,3);
        Ns = -1; % RGBA
    else
        attr = bitset(attr,4);
        Ns = size(M.cdata,2);
    end
end

fwrite(fid,23117,'uint16'); % 0x4D5A, 23117, [77 90], 'MZ'
fwrite(fid,attr,'uint16');
fwrite(fid,Nf,'uint32');
fwrite(fid,Nv,'uint32');
fwrite(fid,0,'uint32');

%-Write data
%--------------------------------------------------------------------------
if Nf > 0, fwrite(fid,M.faces'-1,'uint32');  end
if Nv > 0, fwrite(fid,M.vertices','single'); end
if Ns < 0, fwrite(fid,M.cdata','uint8');     end
if Ns > 0, fwrite(fid,M.cdata,'single');     end

%-Close file
%--------------------------------------------------------------------------
fclose(fid);

%-Optional compression
%--------------------------------------------------------------------------
if fmt
    gzip(filename);
    movefile([filename '.gz'], filename);
end
