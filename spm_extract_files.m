function spm_extract_files(P,cwd)
% FORMAT spm_extract_files(R,cwd)
% forints files (and their subroutines) and expect them to the current
% directory
%__________________________________________________________________________
% Copyright (C) 1994-2011 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_svd.m 4604 2011-12-20 18:21:03Z guillaume $


if nargin == 1; cwd = pwd; end

% deal with cell arrays
%--------------------------------------------------------------------------
if iscell(P)
    for i = 1:length(P)
        spm_extract_files(P{i},cwd)
    end
    return
end

% get file
%--------------------------------------------------------------------------
if isempty(dir(P))
    try
        % check for subroutines
        %------------------------------------------------------------------
        copyfile(which(P),cwd);
    end
else
    return
end


% check for subroutines
%--------------------------------------------------------------------------
try
    fid   = fopen(P);
    Q     = textscan(fid,'%s');
    fclose(fid);
    Q     = Q{1};
    for i = 1:length(Q)
        
        s = strfind(Q{i},'spm_');
        if ~isempty(s)
            j = strfind(Q{i},'(');
            j = j(j > s);
            if ~isempty(j)
                q = [Q{i}(s:(j(1) - 1)) '.m'];
                if ~strcmp(P,q)
                    spm_extract_files(q,cwd);
                end
            end
        end
    end
end

