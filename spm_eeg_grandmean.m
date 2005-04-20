function Do = spm_eeg_grandmean(S)
% average over multiple data sets
% FORMAT Do = spm_eeg_grandmean(S)
% 
% S		    - struct (optional)
% (optional) fields of S:
% P			- filenames (char matrix) of EEG mat-file containing epoched
%             data  
% 
% Output:
% Do		- EEG data struct, result files are saved in the same
%                 directory as first input file.
%_______________________________________________________________________
% 
% spm_eeg_grandmean averages data over multiple ERPs. The data must have
% the same dimensionality and sampling rate. This function can be used for
% grand mean averaging, i.e. computing the average over multiple subjects.
% The output is written to a new file that has the same name as the first
% selected input file, but is prefixed with a 'g'.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% $Id$

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG grandmean setup', 0);

try
    P = S.P;
catch
    P = spm_select(inf, '\.mat$', 'Select EEG mat files');
end

clear D
for i = 1:size(P, 1)
    try
        D{i} = spm_eeg_ldata(deblank(P(i, :)));
    catch    
        error(sprintf('Trouble reading files %s', deblank(P(i, :))));
    end
end

% check dimensionality of data
try
    dim = []; s = [];
    for i = 1:length(D)
        dim(i,:) = size(D{i}.data);
        s(i) = D{i}.Radc;
    end
catch
    error('Data doesn''t have the same dimensionality');
end

if ~all(repmat(dim(1,:), size(dim, 1), 1) == dim)
    error('Data doesn''t have the same dimensionality');
end
if ~all(s(1) == s)
    error('Data doesn''t have the same sampling rate');
end

% test whether data are ERPs
% if D{1}.Nevents ~= D{1}.events.Ntypes
%     error('Data must be ERPs');
% end

% output
Do = D{1};
Do.fnamedat = ['g' D{1}.fnamedat];

fpd = fopen(fullfile(Do.path, Do.fnamedat), 'w');

spm('Pointer', 'Watch'); drawnow;

Do.scale.dim = [1 3];
Do.scale.values = zeros(D{1}.Nchannels, D{1}.Nevents);
	
for i = 1:D{1}.Nevents
        
    d = zeros(D{1}.Nchannels, D{1}.Nsamples);

    for j = 1:D{1}.Nchannels			
        for k = 1:length(D)
            d(j, :) = d(j, :) + D{k}.data(j, :, i);
        end
        d(j, :) = d(j, :)./length(D);
    end
    
    Do.scale.values(:, i) = spm_eeg_write(fpd, d, 2, Do.datatype);

end

fclose(fpd);

Do.data = [];

Do.fname = ['g' D{1}.fname];

D = Do;

if str2num(version('-release'))>=14
    save(fullfile(D.path, D.fname), '-V6', 'D');
else
    save(fullfile(D.path, D.fname), 'D');
end


spm('Pointer', 'Arrow');
