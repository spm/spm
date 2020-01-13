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
% Copyright (C) 2010-2020 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_update.m 7770 2020-01-13 11:26:14Z guillaume $


vspm = spm('Ver');
url  = ['http://www.fil.ion.ucl.ac.uk/spm/download/' lower(vspm) '_updates/'];

if ~nargin
    update = false;
else
    update = true;
end

%-Get list of updates from SPM server
%--------------------------------------------------------------------------
[s,status] = urlread(url);
if ~status
    sts = NaN;
    msg = 'Cannot access SPM server.';
    if ~nargout, error(msg); else varargout = {sts, msg}; end
    return
end

%-Get revision number of latest update
%--------------------------------------------------------------------------
n = regexp(s,[lower(vspm) '_updates_r(\d.*?)\.zip'],'tokens','once');
if isempty(n)
    sts = Inf;
    msg = 'There are no updates available yet.';
    if ~nargout, fprintf([blanks(9) msg '\n']);
    else varargout = {sts, msg}; end
    return
end
n = str2double(n{1});

%-Get revision number of SPM installation
%--------------------------------------------------------------------------
try
    [v,r] = spm('Ver','',1); r = str2double(r);
catch
    error('SPM cannot be found in MATLAB path.');
end
if ~strcmp(v,vspm), error('Your SPM version is %s and not %s',v,vspm); end
rs = [6225 6470 6685 6906 7219 7487 7771];
if isnan(r), r = rs(1); end 
if floor(r) == str2double(vspm(4:end))
    try
        r = rs(round((r-floor(r))*10)+1);
    catch
        r = rs(end);
    end
end

%-Compare versions
%--------------------------------------------------------------------------
if n > r
    sts = n;
    msg = sprintf('         A new version of %s is available on:\n',vspm);
    msg = [msg sprintf('   %s\n',url)];
    msg = [msg sprintf('        (Your version: %d - New version: %d)\n',r,n)];
    if ~nargout, fprintf(msg); else varargout = {sts, msg}; end
else
    sts = 0;
    msg1 = sprintf('Your version of %s is up to date.',vspm);
    msg2 = sprintf('(Your version: %d - Online version: %d)',r,n);
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
        s = unzip([url sprintf('%s_updates_r%d.zip',lower(vspm),n)], d);
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
