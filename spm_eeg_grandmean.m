function Dout = spm_eeg_grandmean(S)
% average over different EEG/MEG data sets
% FORMAT Dout = spm_eeg_average(S)
% 
% S		    - struct (optional)
% (optional) fields of S:
% D			- filenames (char matrix) of EEG mat-file containing epoched
%             data  
% 
% Output:
% Dout			- EEG data struct, result files are saved in the same
%                 directory as first input file.
%_______________________________________________________________________
% 
% spm_eeg_grandmean low-pass filters EEG/MEG epoched data.
%_______________________________________________________________________
% @(#)spm_eeg_grandmean.m	1.1 Stefan Kiebel 04/06/28

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG averaging setup',0);

try
    E = S.D;
catch
    E = spm_get(inf, '.mat', 'Select EEG mat files');
end

P = spm_str_manip(deblank(E(1, :)), 'H');

clear D
for i = 1:size(E, 1)
    try
        D{i} = spm_eeg_ldata(deblank(E(i, :)));
    catch    
        error(sprintf('Trouble reading files %s', deblank(E(i, :))));
    end
end

% output
Dout = D{1};
Dout.fnamedat = ['g' D{1}.fnamedat];

fpd = fopen(fullfile(P, Dout.fnamedat), 'w');

spm('Pointer', 'Watch');

Dout.scale.dim = [1 3];
Dout.scale.values = zeros(D{1}.Nchannels, D{1}.Nevents);
	
for i = 1:D{1}.events.Ntypes
        
    d = zeros(D{1}.Nchannels, D{1}.Nsamples);

    for j = 1:D{1}.Nchannels			
        for k = 1:length(D)
            d(j, :) = d(j, :) + D{k}.data(j, :, i);
        end
        d(j, :) = d(j, :)./length(D);
    end
    
    Dout.scale.values(:, i) = max(abs(d'))./32767;
    d = int16(d./repmat(Dout.scale.values(:, i), 1, Dout.Nsamples));
    fwrite(fpd, d, 'int16');	
end

fclose(fpd);


Dout.data = [];

Dout.fname = ['g' D{1}.fname];
Dout.datatype = 'int16';

D = Dout;
if str2num(version('-release'))>=14
    save(fullfile(P, Dout.fname), '-V6', 'D');
else
    save(fullfile(P, Dout.fname), 'D');
end

spm('Pointer', 'Arrow');
