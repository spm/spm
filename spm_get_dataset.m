function spm_get_dataset(repo, name, rev, outdir)
% Download a dataset from an online repository
% FORMAT spm_get_dataset(repo, name, outdir)
% repo   - name of repository, one of ['spm', 'openfmri']
% name   - name of dataset, e.g. 'auditory' or 'ds000117'
% rev    - revision of dataset [default: '']
% outdir - output directory [default: pwd]
%__________________________________________________________________________
% Copyright (C) 2017-2018 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_get_dataset.m 7265 2018-02-22 15:34:59Z guillaume $


SVNrev = '$Rev: 7265 $';

spm('FnBanner', mfilename, SVNrev);

if nargin < 1, error('A repository name is mandatory.'); end
if nargin < 2, error('A dataset name is mandatory.');    end
if nargin < 3, rev = '';     end
if nargin < 4, outdir = pwd; end

%-Options
%--------------------------------------------------------------------------
timeout = [];
delraw  = true;

%-Get download URL from data repository
%--------------------------------------------------------------------------
switch lower(repo)
    case 'spm'
        base = 'http://www.fil.ion.ucl.ac.uk/spm/download/data/';
        switch lower(name)
            case 'auditory'
                url = fullfile(base,'MoAEpilot','MoAEpilot.zip');
            case 'attention'
                url = fullfile(base,'attention','attention.zip');
            case 'face_rep'
                url = fullfile(base,'face_rep','face_rep.zip');
            case 'face_rfx'
                url = fullfile(base,'face_rfx','face_rfx.zip');
            case 'eeg_mmn'
                url{1} = fullfile(base,'eeg_mmn','subject1.bdf');
                url{2} = fullfile(base,'eeg_mmn','sensors.pol');
            case 'spdcm'
                url = fullfile(base,'spDCM','spDCM.zip');
            otherwise
                error('Unknown dataset "%s" in repository "%s".',name,repo);
        end
        url = strrep(url,'\','/');
        
    case 'openfmri'
        % see https://www.mathworks.com/matlabcentral/answers/92506
        url = 'https://openfmri.org/dataset/api/?format=json';
        [js, sts] = urlread(url,'Timeout',timeout);
        if ~sts, error('Connection to openfMRI failed.'); end
        of = spm_jsonread(js);
        idx = find(ismember({of.accession_number},name));
        if isempty(idx)
            error('Unknown dataset "%s" in repository "%s".',name,repo);
        end
        revs = {of(idx).revision_set.revision_number};
        if isempty(rev)
            rev = revs{end}; % use last revision, sort(revs)?
        else
            if ~ismember(rev,revs)
                warning('Files with revision "%s" not found.',rev);
            end
        end
        lnk = of(idx).link_set;
        url = {};
        for i=1:numel(lnk)
            if strcmp(lnk(i).revision,rev)
                url{end+1} = lnk(i).url;
            end
        end
        
    case 'neurovault'
        % http://neurovault.org/api/?format=api
        error('Work in progress.');
        
    otherwise
        error('Unknown repository "%s".',repo);
end

%-Download
%--------------------------------------------------------------------------
url = cellstr(url);
F   = cell(numel(url),1);
for i=1:numel(url)
    s = get_file_size(url{i});
    f = spm_file(url{i},'filename');
    if ~isnan(s)
        f = [f, sprintf(' [%dM]',round(s/1024/1024))];
    end
    if exist(fullfile(outdir,spm_file(url{i},'filename')),'file') == 7
        F{i} = fullfile(outdir,spm_file(url{i},'filename'));
        sts = true;
        fprintf('%-40s: %30s', f,'...cached');                          %-#
    else
        fprintf('%-40s: %30s', f,'...downloading');                     %-#
        [F{i}, sts] = urlwrite(url{i}, ...
            fullfile(outdir,spm_file(url{i},'filename')), 'Timeout',timeout);
    end
    if sts, msg = '...done'; else, msg = '...failed';end
    fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),msg);                 %-#
    if ~sts, error('Download failed.'); end
end
    
%-Uncompress (and delete archive)
%--------------------------------------------------------------------------
filenames = {};
for i=1:numel(F)
    switch lower(spm_file(F{i},'ext'))
        case 'zip'
            f = unzip(F{i},spm_file(F{i},'path'));
            if delraw, spm_unlink(F{i}); end
            filenames = [filenames; f(:)];
        case {'tar','gz','tgz'}
            f = untar(F{i},spm_file(F{i},'path'));
            if delraw, spm_unlink(F{i}); end
            filenames = [filenames; f(:)];
        otherwise
            filenames = F;
    end
end

fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#


%==========================================================================
% function s = get_file_size(url)
%==========================================================================
function s = get_file_size(url)
try
    uri    =  matlab.net.URI(url);
    method = matlab.net.http.RequestMethod.HEAD;
    req    = matlab.net.http.RequestMessage(method);
    resp   = req.send(uri);
    s      = convert(getFields(resp,'Content-Length'));
    if isempty(s), s = NaN; end
catch
    s = NaN;
end
