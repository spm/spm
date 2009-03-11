function D = spm_eeg_robust_average_TF(S)
% Robust averaging TF data based on weights calculated on ERP of same data.
% FORMAT D = spm_eeg_robust_average_TF(S)
%
% S        - optional input struct
% (optional) fields of S:
%   S.D    - 
%   S.D2   - 
%   S.c    - 
%
% D        - MEEG object (also written to disk)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% James Kilner
% $Id: spm_eeg_robust_average_TF.m 2861 2009-03-11 18:41:03Z guillaume $

SVNrev = '$Rev: 2861 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','M/EEG robust average TF'); spm('Pointer','Watch');

%-Get MEEG objects
%--------------------------------------------------------------------------
try
    D = S.D;
catch
    [D, sts] = spm_select(1, 'mat', 'Select M/EEG mat file');
    if ~sts, D = []; return; end
    S.D = D;
end

D = spm_eeg_load(D);

try
    D2 = S.D2;
catch
    [D2, sts] = spm_select(1, 'mat', 'Select robust averaged erp file file');
    if ~sts, return; end
    S.D2 = D2;
end

D2 = spm_eeg_load(D2);
P  = spm_str_manip(D2, 'H');

%-Get parameters
%--------------------------------------------------------------------------
try
    c = S.c;
catch
    % if there is no S.c, assume that user wants default average within
    % trial type
    c = eye(D.events.Ntypes);
end

%-
%--------------------------------------------------------------------------
d=zeros(size(D.data,1),size(D.data,2),size(D.data,3),D.events.Ntypes);
fh=fopen(fullfile(D.path,D.fnamedat),'r');
  spm_progress_bar('Init', D.Nevents, 'averaging'); drawnow;
    if D.Nevents > 100, Ibar = floor(linspace(1, D.Nevents,100));
    else, Ibar = [1:D.Nevents]; end
for n=1:D.Nevents
  
    if D.events.reject(n)== 0
        [m,i]=find(D.events.types==D.events.code(n));
        data=fread(fh,[1,size(D.data,1)*size(D.data,2)*size(D.data,3)],'short');
        data=reshape(data,size(D.data,1),size(D.data,2),size(D.data,3),1);
        data=data.*repmat(D.scale(:,1,1,n),[1,D.Nfrequencies, D.Nsamples]);
        for f=1:D.Nfrequencies
            data(:,f,:)=squeeze(data(:,f,:)).*D2.weights(:,(n-1)*size(D.data,3)+1:n*size(D.data,3));
        end
        d(:,:,:,i)=d(:,:,:,i)+data;
    end
    
    if ismember(n, Ibar)
        spm_progress_bar('Set', n);
        drawnow;
    end
    
    
end


spm_progress_bar('Clear');
D.fnamedat = ['ma' D.fnamedat];

fpd = fopen(fullfile(P, D.fnamedat), 'w');
D.scale = zeros(D.Nchannels, 1, 1, D.events.Ntypes);

for n=1:D.events.Ntypes
    dat=squeeze(d(:,:,:,n)./length(find(D.events.code==D.events.types(n) & ~D.events.reject)));
    D.scale(:, 1, 1, n) = max(max(squeeze(abs(dat)), [], 3), [], 2)./32767;
    dat = int16(dat./repmat(D.scale(:, 1, 1, n), [1, D.Nfrequencies, D.Nsamples]));
    fwrite(fpd, dat, 'int16');
end
fclose (fh) ;

fclose(fpd);
D.Nevents = size(c, 2);

% labeling of resulting contrasts, take care to keep numbers of old trial
% types
% check this again: can be problematic, when user mixes within-trialtype
% and over-trial type contrasts
D.events.code = size(1, size(c, 2));
for i = 1:size(c, 2)
    if sum(c(:, i)) == 1 & sum(~c(:, i)) == size(c, 1)-1
        D.events.code(i) = find(c(:, i));
    else
        D.events.code(i) = i;
    end
end

D.events.time = [];
D.events.types = D.events.code;
D.events.Ntypes = length(D.events.types);
D.data = [];
D.events.reject = zeros(1, D.Nevents);
D.events.blinks = zeros(1, D.Nevents);
D.fname = ['ma' D.fname];
D.weights=D2.weights;
if spm_matlab_version_chk('7') >= 0
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

%-Cleanup
%--------------------------------------------------------------------------
spm_progress_bar('Clear');
spm('FigName','M/EEG robust average TF: done'); spm('Pointer','Arrow');
