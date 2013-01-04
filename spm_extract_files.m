function spm_extract_files(P,cwd)
% FORMAT spm_extract_files(P,cwd)
% forints files (and their subroutines) and expect them to the current
% directory
%__________________________________________________________________________
% Copyright (C) 2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_extract_files.m 5175 2013-01-04 12:50:44Z guillaume $


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
        for k = 1:length(s)
            
            % calls: spm_???(
            %--------------------------------------------------------------
            j = strfind(Q{i},'(');
            j = j(find(j > s(k),1));
            if ~isempty(j)
                q = [Q{i}(s(k):(j - 1)) '.m'];
                if ~strcmp(P,q)
                    spm_extract_files(q,cwd);
                end
            end
            
            % functions: 'spm_???'
            %--------------------------------------------------------------
            j = strfind(Q{i},'''');
            j = j(find(j > s(k),1));
            if ~isempty(j)
                q = [Q{i}(s(k):(j - 1)) '.m'];
                if ~strcmp(P,q)
                    spm_extract_files(q,cwd);
                end
            end
            
        end
    end
end
