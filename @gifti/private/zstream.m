function Z = zstream(action,D)
% Compress/decompress stream of bytes using Deflate/Inflate
% FORMAT Z = zstream('C',D)
% D        - data stream to compress (converted to uint8 if needed)
% Z        - compressed data stream (uint8 vector)
% FORMAT D = zstream('D',Z)
% Z        - data stream to decompress (uint8 vector)
% D        - decompressed data stream (uint8 vector)
%
% If action is upper case ('C','D'), a zlib stream is used (zlib header
% with an adler32 checksum). Otherwise, if action is lower case ('c','d'),
% a raw deflate stream is assumed.
%__________________________________________________________________________
%
% This C-MEX file relies on:
% * miniz, by Rich Geldreich
%   https://github.com/richgel999/miniz
% Fallback Java implementation is adapted from:
% * dzip/dunzip, by Michael Kleder
%   https://www.mathworks.com/matlabcentral/fileexchange/8899
%__________________________________________________________________________
% Copyright (C) 2015-2018 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: zstream.m 7399 2018-08-16 12:04:14Z guillaume $


if exist('OCTAVE_VERSION','builtin')
    error('zstream.c not compiled - see Makefile');
end

switch upper(action)
    case 'C'
        D = typecast(D(:),'uint8');
        f = java.io.ByteArrayOutputStream();
        g = java.util.zip.DeflaterOutputStream(f,java.util.zip.Inflater(action~='C'));
        g.write(D);
        g.close;
        Z = typecast(f.toByteArray,'uint8');
        f.close;
        
    case 'D'
        import com.mathworks.mlwidgets.io.InterruptibleStreamCopier
        a   = java.io.ByteArrayInputStream(D);
        b   = java.util.zip.InflaterInputStream(a,java.util.zip.Inflater(action~='D'));
        isc = InterruptibleStreamCopier.getInterruptibleStreamCopier;
        c   = java.io.ByteArrayOutputStream;
        isc.copyStream(b,c);
        Z   = c.toByteArray;
        
    otherwise
        error('Unknown action "%s".',action);
end
