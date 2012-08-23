function spm_update(update)
% Check (and install) SPM updates from the FIL server
% FORMAT spm_update
% This function will connect itself to the FIL server, compare the
% version number of the updates with the one of the SPM installation 
% currently in the MATLAB path and will display the outcome.
%
% FORMAT spm_update(update)
% Invoking this function with any input parameter will do the same as
% above but will also attempt to download and install the updates.
%__________________________________________________________________________
% Copyright (C) 2010-2012 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_update.m 4857 2012-08-23 17:17:30Z guillaume $

url  = 'ftp://ftp.fil.ion.ucl.ac.uk/spm/spm12_updates/';
vspm = 'SPM12';

if ~nargin
    update = false;
else
    update = true;
end

[s,sts] = urlread(url);
if ~sts, error('Cannot access SPM server.'); end
n       = regexp(s,[lower(vspm) '_updates_r(\d.*?)\.zip'],'tokens','once');
if isempty(n)
    fprintf('         There are no updates available yet.\n');
    return;
else
    n   = str2double(n{1});
end

try
    [v,r] = spm('Ver','',1); r = str2double(r);
catch
    error('SPM cannot be found in MATLAB path.');
end
if ~strcmp(v,vspm), error('Your SPM version is %s and not %s',v,vspm); end
rs = [];
if isnan(r), r = rs(1); end 
if floor(r) == str2double(vspm(4:end))
    try
        r = rs(round((r-floor(r))*10)+1);
    catch
        r = rs(end);
    end
end

if n > r
    fprintf('         A new version of %s is available on:\n',vspm);
    fprintf('     %s\n',url);
    fprintf('        (Your version: %d - New version: %d)\n',r,n);

    if update
        d = spm('Dir'); 
        delete(get(0,'Children')); spm('clean'); evalc('spm_rmpath');
        try
            lastwarn('');
            s = unzip([url sprintf('%s_updates_r%d.zip',lower(vspm),n)], d);
            fprintf('             %d files have been updated.\n',numel(s));
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
else
    fprintf('         Your version of %s is up to date.\n',vspm);
end
