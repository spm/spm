function V = spm_filename_check(V)
% Checks paths are valid and tries to restore path names
% FORMAT V = spm_filename_check(V)
%
% V - cell array for file handles
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_filename_check.m 2708 2009-02-06 19:52:48Z karl $

% check filenames
%--------------------------------------------------------------------------
for i = 1:length(V)


    % see if file exists
    %----------------------------------------------------------------------
    if ~spm_existfile(V(i).fname)

        % try current directory
        %------------------------------------------------------------------
        [p,n,e] = fileparts(V(i).fname);
        fname   = which([n,e]);
        if ~isempty(fname)
            V(i).fname = fname;
        else

            % try parent directory
            %--------------------------------------------------------------
            cwd = pwd;
            cd('..')
            fname  = which([n,e]);
            if ~isempty(fname)
                V(i).fname = fname;
            else

                % try children of parent
                %----------------------------------------------------------
                V = spm_which_filename(V);
                cd(cwd)
                return

            end
            cd(cwd)
        end

    end

end


function V = spm_which_filename(V)


% get children directories of parent
%--------------------------------------------------------------------------
cwd = pwd;
cd('..')
gwd = genpath(pwd);
addpath(gwd);

% cycle through handles
%--------------------------------------------------------------------------
for i = length(V)

    try
        % get relative path (directory and filename) and find in children
        %------------------------------------------------------------------
        j     = findstr(V(i).fname,filesep);
        fname = which(fname(j(end - 1):end));
        if ~isempty(fname)
            V(i).fname = fname;
        end
    end
end

% reset paths
%--------------------------------------------------------------------------
rmpath(gwd);
cd(cwd);

