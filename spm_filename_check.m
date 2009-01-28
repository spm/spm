function V = spm_filename_check(V)
% Checks paths are valid and tries to restore path names
% FORMAT V = spm_filename_check(V)
%
% V - cell array for file handles
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_filename_check.m 2664 2009-01-28 20:25:20Z karl $

% check filenames
%--------------------------------------------------------------------------
for i = 1:length(V)

    % try current directory
    %----------------------------------------------------------------------
    [p,n,e] = fileparts(V(i).fname);
    fname   = [n,e];
    try
        V(i).fname = which(fname);
    catch

        % try parent directory
        %------------------------------------------------------------------
        cwd = pwd;
        cd('..')
        try
            V(i).fname = which(fname);
        catch

            % try children of parent
            %--------------------------------------------------------------
            V = spm_which_filename(V);
            return

        end
        cd(cwd)
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
        j          = findstr(V(i).fname,filesep);
        V(i).fname = which(fname(j(end - 1):end));
    end

end

% reset paths
%--------------------------------------------------------------------------
rmpath(gwd);
cd(cwd);

