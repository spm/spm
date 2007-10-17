function D = spm_eeg_weight_epochs(S);
% computes contrasts over trials or trial types.
% FORMAT D = spm_eeg_weight_epochs(S)
%
% S		    - optional input struct
% (optional) fields of S:
% D			- filename of EEG mat-file with epoched data
% c         - contrast matrix, each row computes a contrast of the data
%
% Output:
% D			- EEG data struct (also written to files)
%_______________________________________________________________________
%
% spm_eeg_weight_epochs computes contrasts of data, over epochs of data. The
% input is a single MEEG file.
% The argument c must have dimensions N_contrasts X N_epochs, where N_contrasts is
% the number of contrasts and N_epochs the number of epochs, i.e. each row of c
% contains one contrast vector. The output
% is a MEEG file with N_contrasts epochs. The typical use is to compute,
% for display purposes, contrasts like the difference or interaction
% between trial types in channel space.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Stefan Kiebel
% Extended by Rik Henson to handle weighted averages and TF
% $Id: spm_eeg_weight_epochs.m 955 2007-10-17 15:15:09Z rik $

[Finter,Fgraph,CmdLine] = spm('FnUIsetup','EEG averaging setup',0);

try
    D = S.D;
catch
    D = spm_select(1, '.*\.mat$', 'Select EEG mat file');
end

P = spm_str_manip(D, 'H');

try
    D = spm_eeg_ldata(D);
catch
    error(sprintf('Trouble reading file %s', D));
end

try
    c = S.c;
catch
    c = spm_input('Enter contrasts', '1', 'x', '', inf, eye(D.Nevents));
end

if isfield(D.events,'repl')
  try
    WeightAve = S.WeightAve;
  catch
    WeightAve = spm_input('Weight by num replications?', '+1', 'yes|no', [1 0]);
  end
else
  WeightAve = 0;
end

% incorporate sanity check of c, later

D.fnamedat = ['m' D.fnamedat];

fpd = fopen(fullfile(P, D.fnamedat), 'w');

spm('Pointer', 'Watch'); drawnow;

N_contrasts = size(c, 1);

if isfield(D, 'Nfrequencies')
  D.scale = zeros(D.Nchannels, 1, 1, N_contrasts);
  D.datatype = 'int16';
else
  D.scale = zeros(D.Nchannels, 1, N_contrasts);
end

spm_progress_bar('Init', N_contrasts, 'Contrasts computed'); drawnow;
if N_contrasts > 100, Ibar = floor(linspace(1, N_contrasts, 100));
else, Ibar = [1:N_contrasts]; end

for i = 1:N_contrasts
    
  if WeightAve
    	p = find(c(i,:)==1);
    	if ~isempty(p)
        	r = D.events.repl(p);
        	c(i,p) = r/sum(r);
    	end

   	p = find(c(i,:)==-1);
    	if ~isempty(p)
        	r = D.events.repl(p);
        	c(i,p) = -r/sum(r);
    	end
  end
  disp(['Contrast ',mat2str(i),': ',mat2str(c(i,:),3)])

  if isfield(D, 'Nfrequencies')
	
    d = zeros(D.Nchannels, D.Nfrequencies, D.Nsamples);

    % Could be made faster!!!
    for j = 1:D.Nchannels
      for f = 1:D.Nfrequencies
        d(j, f, :) = c(i,:)*squeeze(D.data(j, f, :, :))';
      end
    end

    D.scale(:, 1, 1, i) = max(max(abs(d), [], 3), [], 2)./32767;
    d = int16(d./repmat(D.scale(:, 1, 1, i), [1, D.Nfrequencies, D.Nsamples]));

%    D.scale(:, 1, 1, i) = (max(abs(reshape(d, [D.Nchannels D.Nfrequencies*D.Nsamples])'))./32767);
%    d = int16(round(d./repmat(D.scale(:, 1, 1, i), [1 D.Nfrequencies D.Nsamples])));
    fwrite(fpd, d, 'int16');	

  else

    d = zeros(D.Nchannels, D.Nsamples);

    %if D.events.repl(i) == 0
    %  warning('%s: No trials for contrast %d', D.fname, i);
    % else

    for j = 1:D.Nchannels
        d(j, :) = c(i,:)*squeeze(D.data(j, :, :))';
    end
    %end

    D.scale(:, 1, i) = spm_eeg_write(fpd, d, 2, D.datatype);

  end

  if ismember(i, Ibar)
        spm_progress_bar('Set', i);
        drawnow;
  end

  newrepl(i) = sum(D.events.repl(find(c(i,:)~=0)));

end

spm_progress_bar('Clear');

fclose(fpd);

D.Nevents = N_contrasts;
D.events.code = [1:N_contrasts];
D.events.repl = newrepl;
D.events.contrast = c;

D.events.time = [];
D.events.types = D.events.code;
D.events.Ntypes = N_contrasts;

D.data = [];
D.events.reject = zeros(1, D.Nevents);
D.events.blinks = zeros(1, D.Nevents);

D.fname = ['m' D.fname];

if spm_matlab_version_chk('7') >= 0,
    save(fullfile(P, D.fname), '-V6', 'D');
else
    save(fullfile(P, D.fname), 'D');
end

spm('Pointer', 'Arrow');
