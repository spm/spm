function bool = bf_isfield(BF,field)
% Efficiently identify if a field is contained within a BF file
% FORMAT bool = bf_isfield(BF,field)
%__________________________________________________________________________

% George O'Neill
% Copyright (C) 2020-2023 Wellcome Centre for Human Neuroimaging

if exist(BF,'file') || exist(BF,'var')
    if isstruct(BF)
        bool = isfield(BF,field);
    elseif ischar(BF)
        bool = false;
        fid = H5F.open(BF); % Use low level H5F builtin to open
        try
            gid = H5G.open(fid,['/' field]);
            hInfo = H5G.get_info(gid);
            bool = hInfo.nlinks > 0;
            H5G.close(gid);
        end
        H5F.close(fid);
    else
        error('BF must be either a struct or string.')
    end
else
    bool = false;
end
