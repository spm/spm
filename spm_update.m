function varargout = spm_update(update)
% Check (and install) SPM updates from the FIL server
% FORMAT spm_update
% This function will contact the FIL server, compare the version number of
% the updates with the one of the SPM installation currently in the MATLAB
% path and display the outcome.
%
% FORMAT spm_update(update)
% Invoking this function with any input parameter will do the same as
% above but will also attempt to download and install the updates.
% Note that it will close any open window and clear the workspace.
%
% FORMAT [sts, msg] = spm_update(update)
% sts  - status code:
%        NaN - SPM server not accessible
%        Inf - no updates available
%        0   - SPM installation up to date
%        n   - new revision <n> is available for download
% msg  - string describing outcome, that would otherwise be displayed.
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2010-2022 Wellcome Centre for Human Neuroimaging


[spm_version, spm_revision]  = spm('Ver');

if ~nargin
    update = false;
else
    update = true;
end

%-Get list of updates from SPM server
%--------------------------------------------------------------------------
try
    response = webread('https://api.github.com/repos/spm/spm/releases');
catch
    sts = NaN;
    msg = 'Cannot access SPM server.';
    if ~nargout, error(msg); else varargout = {sts, msg}; end
    return
end

%-Get latest version
%--------------------------------------------------------------------------
valid_version_pattern = '^\d{2}\.\d{2}(\.\d+)?$';

if iscell(response)
    tagged_versions = string(cellfun(@(r) r.tag_name, response, 'uni', 0));
else 
    tagged_versions = string({response.tag_name});
end

valid_versions = ~cellfun('isempty', regexp(tagged_versions, valid_version_pattern));
sorted_versions = sort(tagged_versions(valid_versions), 'descend');

if numel(sorted_versions) == 0
    sts = Inf;
    msg = 'There are no updates available yet.';
    if ~nargout, fprintf([blanks(9) msg '\n']);
    else varargout = {sts, msg}; end
    return
end

latest_version = sorted_versions{1};
url = sprintf('https://github.com/spm/spm/releases/download/%s/spm_%s.zip', latest_version, latest_version);

%-Compare versions
%--------------------------------------------------------------------------
if string(latest_version) > string(spm_revision)
    sts = spm_revision;
    msg = sprintf('         A new version of SPM is available on:\n');
    msg = [msg sprintf('   %s\n',url)];
    msg = [msg sprintf('        (Your version: %s - New version: %s)\n',spm_revision,latest_version)];
    if ~nargout, fprintf(msg); else varargout = {sts, msg}; end
else
    sts = 0;
    msg1 = sprintf('Your version of %s is up to date.',spm_version);
    msg2 = sprintf('(Your version: %s - Online version: %s)',spm_revision,latest_version);
    if ~nargout, fprintf([blanks(9) msg1 '\n' blanks(5) msg2 '\n']);
    else varargout = {sts, sprintf('%s\n%s',msg1,msg2)}; end
    return
end

%-and install...
%--------------------------------------------------------------------------
if update
    d = spm('Dir');
    delete(get(0,'Children')); spm('clean'); evalc('spm_rmpath'); drawnow
    try
        lastwarn('');
        m = '          Download and install in progress...\n';
        if ~nargout, fprintf(m); else varargout = {sts, [msg m]}; end
        s = unzip(url, [d '/spm_update_tmp']);
        update_folder = [d '/spm_update_tmp/spm'];
        movefile(fullfile(update_folder, '*'), d);
        rmdir(update_folder, 's');
        m = sprintf('         Success: %d files have been updated.\n',numel(s));
        if ~nargout, fprintf(m); else varargout = {sts, [msg m]}; end
    catch
        le = lasterror;
        switch le.identifier
            case 'MATLAB:checkfilename:urlwriteError'
                fprintf('          Update failed: cannot download update file.\n');
            otherwise
                fprintf('\n%s\n',le.message);
        end
    end
    [warnmsg, msgid] = lastwarn;
    switch msgid
        case ''
        case 'MATLAB:extractArchive:unableToCreate'
            fprintf('          Update failed: check folder permission.\n');
        case 'MATLAB:extractArchive:unableToOverwrite'
            fprintf('          Update failed: check file permissions.\n');
        otherwise
            fprintf('          Update failed: %s.\n',warnmsg);
    end
    addpath(d);
    rehash toolboxcache;
end
