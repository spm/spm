function spm_get_dataset(store, name, outdir)
% Download datasets
% FORMAT spm_get_dataset(store, name, outdir)
% store  - one of ['spm', 'openfmri']
% name   - name of dataset, e.g. 'auditory' or 'ds000117'
% outdir - output directory [default: pwd]
%__________________________________________________________________________
% Copyright (C) 2017 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_get_dataset.m 7138 2017-07-21 10:57:37Z guillaume $


SVNrev = '$Rev: 7138 $';

spm('FnBanner', mfilename, SVNrev);

if nargin < 3
    outdir = pwd;
end

timeout = [];

%-Get download URL from datastore
%--------------------------------------------------------------------------
switch lower(store)
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
            otherwise
                error('Unknown dataset "%s" in datastore "%s".',name,store);
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
            error('Unknown dataset "%s" in datastore "%s".',name,store);
        end
        revs = {of(idx).revision_set.revision_number};
        rev = revs{end}; % choose revision
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
        error('Unknown datastore "%s".',store);
end

%-Download
%--------------------------------------------------------------------------
url = cellstr(url);
F   = cell(numel(url),1);
for i=1:numel(url)
    [F{i}, sts] = urlwrite(url{i}, ...
        fullfile(outdir,spm_file(url{i},'filename')), 'Timeout',timeout);
    if ~sts, error('Download failed.'); end
end
    
%-Uncompress (and delete archive)
%--------------------------------------------------------------------------
filenames = {};
for i=1:numel(F)
    switch lower(spm_file(F{i},'ext'))
        case 'zip'
            f = unzip(F{i},spm_file(F{i},'path'));
            spm_unlink(F{i});
            filenames = [filenames; f(:)];
        case {'tar','gz','tgz'}
            f = untar(F{i},spm_file(F{i},'path'));
            spm_unlink(F{i});
            filenames = [filenames; f(:)];
        otherwise
            filenames = F;
    end
end

fprintf('%-40s: %30s\n','Completed',spm('time'));                       %-#
