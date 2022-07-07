function V = spm_check_filename(V)
% Check paths are valid and try to restore path names
% FORMAT V = spm_check_filename(V)
%
% V - struct array of file handles
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2009-2022 Wellcome Centre for Human Neuroimaging

if isdeployed, return; end

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


%==========================================================================
function V = spm_which_filename(V)


% get children directories of parent
%--------------------------------------------------------------------------
cwd = pwd;
cd('..')
gwd = genpath(pwd);
addpath(gwd);

% cycle through handles
%--------------------------------------------------------------------------
for i = 1:length(V)

    try
        % get relative path (directory and filename) and find in children
        %------------------------------------------------------------------
        j     = strfind(V(i).fname,filesep);
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
