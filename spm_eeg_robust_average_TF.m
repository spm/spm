function D = spm_eeg_robust_average_TF(S)
% robust averaging TF data based on weights calculated on ERP of same data.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%James Kilner
% $Id: spm_eeg_artefact.m 265 2005-10-19 17:24:54Z james $


[Finter,Fgraph,CmdLine] = spm('FnUIsetup', 'EEG artefact setup',0);

try
    D = S.D;
catch
    D = spm_select(1, '.*\.mat$', 'Select TF mat file');
end

P = spm_str_manip(D, 'H');

try
    D = spm_eeg_ldata(D);
catch
    error(sprintf('Trouble reading file %s', D));
end
try
    D2 = S.D2;
catch
    D2 = spm_select(1, '.*\.mat$', 'Select robust averaged erp file file');
end

P = spm_str_manip(D2, 'H');

try
    D2 = spm_eeg_ldata(D2);
catch
    error(sprintf('Trouble reading file %s', D));
end
try
    c = S.c;
catch
    % if there is no S.c, assume that user wants default average within
    % trial type
    c = eye(D.events.Ntypes);
end


spm('Clear',Finter, Fgraph);

[Finter,Fgraph,CmdLine] = spm('FnUIsetup', 'EEG artefact setup',0);


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
fclose (fh)	;

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