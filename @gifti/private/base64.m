function y = base64(action,x)
% Base64 binary-to-text encoding/decoding scheme
% FORMAT y = zstream('encode',x)
% x        - data stream to encode (uint8)
% y        - Base64-encoded data stream (uint8)
% FORMAT y = zstream('decode',x)
% x        - data stream to decode (uint8)
% y        - Base-64 decoded data stream (uint8)
%__________________________________________________________________________
%
% This C-MEX file is a wrapper around:
%   https://stackoverflow.com/a/37109258
% by polfosol: https://stackoverflow.com/users/5358284/polfosol
%
% >> char(base64('decode',base64('encode',uint8('Base64'))))
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2008-2023 Wellcome Centre for Human Neuroimaging


switch action
    case 'encode'
        y = base64encode(x);
    case 'decode'
        y = base64decode(x);
    otherwise
        error('Unknown action.');
end
