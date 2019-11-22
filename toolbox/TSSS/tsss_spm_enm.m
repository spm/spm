function D = tsss_spm_enm(S)
% Perform tSSS on rawdata file 'infile' acquired with an Elekta Neuromag
% 306-channel MEG system. The tSSS-processed data are written into 'outfile'.
% The SSS operation is performed in the head coordinate system with the
% expansion origin given by the 3x1 dimensional vector 'r_sphere' ([x y z]
% m; typically [0 0 0.04]). The temporal correlation analysis of tSSS is based
% on raw data segments of length 't_window' (in seconds) and correlation limit
% 'corr_limit'. The order values of the internal and external SSS bases are
% 'Lin' and 'Lout', typically 8 and 3, respectively. 
%
% NOTE: This tSSS function does not utilize the so-called fine calibration
% information of the MEG system. Also, no basis vector selection is
% performed as a regularization step. Ideally, the input file 'infile'
% should contain cross-talk compensated data, which can be done by the
% MaxFilter software.

% Copyright (c) 2016, Elekta Oy
% ---------------------------------------
% 
% Redistribution and use of the Software in source and binary forms, with or without 
% modification, are permitted for non-commercial use.
% 
% The Software is provided "as is" without warranties of any kind, either express or
% implied including, without limitation, warranties that the Software is free of defects,
% merchantable, fit for a particular purpose. Developer/user agrees to bear the entire risk 
% in connection with its use and distribution of any and all parts of the Software under this license.
% 

% $Id: tsss_spm_enm.m 7703 2019-11-22 12:06:29Z guillaume $

SVNrev = '$Rev: 7703 $';

%-Startup
%--------------------------------------------------------------------------
spm('FnBanner', mfilename, SVNrev);
spm('FigName','TSSS tool'); spm('Pointer','Watch');

if ~isfield(S, 'tsss'),            S.tsss       = 1;     end
if ~isfield(S, 'refind'),          S.refind     = 1;     end
if ~isfield(S, 'Dref'),            S.Dref       = '';    end
if ~isfield(S, 't_window'),        S.t_window   = 1;     end
if ~isfield(S, 'corr_limit'),      S.corr_limit = 0.98;  end
if ~isfield(S, 'magscale'),        S.magscale   = 100;   end
if ~isfield(S, 'xspace'),          S.xspace     = 0;     end
if ~isfield(S, 'Lin'),             S.Lin        = 8;     end
if ~isfield(S, 'Lout'),            S.Lout       = 3;     end
if ~isfield(S, 'cond_threshold'),  S.cond_threshold = 50; end
if ~isfield(S, 'prefix'),          S.prefix     = 'sss_'; end

D = spm_eeg_load(S.D);

isneuromag = strncmpi(ft_senstype(D.chanlabels), 'neuromag', 7);

if ~isempty(S.Dref)
    Dref      = spm_eeg_load(S.Dref);
    Dref      = montage(Dref, 'switch', 0);
    realign   = true;    
    ref.grad  = sensors(Dref, 'MEG');
    ref.fid   = fiducials(Dref);
    
    if isfield(Dref, 'SSS')
        ref.SSS = Dref.SSS(S.refind);
    else
        error('Run SSS on the reference dataset first');
    end
    
    meg_ch = Dref.indchantype({'MEG', 'MEGPLANAR'});
    if isneuromag
        mags    = strmatch('MEGMAG', Dref.chantype(meg_ch));
        grads   = strmatch('MEGPLANAR', Dref.chantype(meg_ch));
    else
        mags     = 1:length(meg_ch);
        grads    = [];
    end
    
    ref.meg_ch = meg_ch;
    ref.mags    = mags;
    ref.grads   = grads;
else
    realign   = S.refind>0;
    ref = [];
    Dref = D;
end

mag2SI = 1e-15;
grad2SI = 1e-12;

tsss = S.tsss;

r_sphere   = [0 0 0.04]';
t_window   = S.t_window;
corr_limit = S.corr_limit;
Lin  = S.Lin;
Lout = S.Lout;


magscale = S.magscale;%;100;

if isequal(D.type, 'continuous')
    cont = 1;
    nsamp = D.nsamples;
    data_length = nsamp/D.fsample;
    nwin = floor(data_length/t_window);
    nsamp_win = floor(nsamp/nwin);  % Number of samples in one window
else
    cont = 0;
    nsamp_win = D.nsamples;
    nwin = D.ntrials;
end


meg_ch  = D.indchantype({'MEG', 'MEGPLANAR'});
other_ch = setdiff(1:D.nchannels, meg_ch);
if isneuromag
    mags    = strmatch('MEGMAG', D.chantype(meg_ch));
    grads   = strmatch('MEGPLANAR', D.chantype(meg_ch));
else
    mags     = 1:length(meg_ch);
    grads    = [];
end
%%
goodchs = find(~D.badchannels(meg_ch));
%%
coils = zeros(1, length(meg_ch));
coils(mags) = 1;
%
% Calculate SSS basis
%
if isneuromag
    origheader = D.origheader;
    if ~iscell(origheader)
        origheader = {origheader};
    end
    
    nlocations = numel(origheader);
else
    if isfield(D, 'allsens')
        allsens = D.allsens;
        allfid  = D.allfid;
    else
        allsens = D.sensors('MEG');
        allfid  = D.fiducials;
    end
    
    nlocations = numel(allsens);
end

for i = 1:nlocations
        
    disp('Calculating the SSS basis...');
    if isneuromag
        [R,EX,EY,EZ] = origheader_getpos(origheader{i}, 'head');%
        Sin = Sin_vsh_fast(r_sphere,R,EX,EY,EZ,coils,Lin);
    else       
        [R,EX,EY,EZ] = ft_getpos(allsens(i), allfid(i), D.chanlabels(meg_ch));%
        %[R,EX,EY,EZ] = origheader_getpos(D.origheader, 'head');%
         Sin = Sin_vsh_ctf(r_sphere,R,EX,EY,EZ,Lin);
    end
    Sin(mags,:) = magscale*Sin(mags,:);
    if Lout>0
        if isneuromag
            Sout = Sout_vsh_fast(r_sphere,R,EX,EY,EZ,coils,Lout);
        else
            Sout = Sout_vsh_ctf(r_sphere,R,EX,EY,EZ,Lout);
        end
        Sout(mags,:) = magscale*Sout(mags,:);
    else
        Sout  = [];
        SNout = [];
    end    
    %%
    
    for j = 1:size(Sin,2)
        SNin(:,j) = Sin(:,j)/norm(Sin(:,j));
    end    
    for j = 1:size(Sout,2)
        SNout(:,j) = Sout(:,j)/norm(Sout(:,j));
    end
    SN = [SNin SNout];
    pSN = pinv(SN(goodchs,:));
    SSS(i) = struct('Sin', Sin, 'SNin', SNin, 'Sout', Sout, 'SNout', SNout, 'SN', SN, 'pSN', pSN);
end


out_meg_ch = meg_ch;
out_other_ch = other_ch;
out_mags = mags;
out_grads = grads;
if realign
    if ~isempty(ref)
        out_meg_ch = 1:length(ref.meg_ch);
        out_other_ch = length(out_meg_ch)+[1:length(other_ch)];
        out_mags   = ref.mags;
        out_grads  = ref.grads;
    else
        if isfield(D, 'allsens')
            ref.grad  = D.allsens(S.refind);
            ref.fid   = D.allfid(S.refind);            
        else
            ref.grad = D.sensors('MEG');
            ref.fid  = D.fiducials;
        end
        ref.meg_ch = meg_ch;
        ref.mags   = mags;
        ref.grads  = grads;
        ref.SSS    = SSS(S.refind);
    end
end

nchanout = length(out_meg_ch)+length(out_other_ch);

Dout = clone(D, spm_file(fullfile(D), 'prefix', S.prefix), [nchanout D.nsamples D.ntrials], 1);

if isfield(D, 'fileind')
    fileind = D.fileind;
else
    fileind = ones(1, D.ntrials);
end

spm_progress_bar('Init', nwin, 'Windows processed');
if nwin > 100, Ibar = floor(linspace(1, nwin, 100));
else Ibar = [1:nwin]; end

%%
%
% Process the raw data
%

for j = 1:nwin
    spm_progress_bar('Set','ylabel','reading...');
    if cont
        from = (j-1)*nsamp_win + 1;
        if j<nwin
            to = j*nsamp_win;
        else
            to = nsamp;
        end
        data = D(:, from:to);
        
        bad  = badsamples(D, meg_ch, from:to, 1);
    else
        data = D(:, :, j);
        
        bad  = badsamples(D, meg_ch, ':', j);
        
        
        SNin  = SSS(fileind(j)).SNin;
        SNout = SSS(fileind(j)).SNout;
        SN    = SSS(fileind(j)).SN;
        Sin   = SSS(fileind(j)).Sin;
        if ~isequal(pSN, SSS(fileind(j)).pSN)
            pSN   = SSS(fileind(j)).pSN;
            goodchs = find(~D.badchannels(meg_ch));
        end
    end
    
    cgoodchs = find(~any(bad'));
    if ~isequal(cgoodchs, goodchs);
        goodchs = cgoodchs;
        pSN = pinv(SN(goodchs,:));
    end
    
    if realign
        Sin  = ref.SSS.Sin;       
    end
    
    spm_progress_bar('Set','ylabel','processing...');
    data(meg_ch(mags),:) = magscale*mag2SI*data(meg_ch(mags),:);
    data(meg_ch(grads),:) = grad2SI*data(meg_ch(grads),:);
    if size(data,2) >= nsamp_win
        if norm(data(meg_ch(goodchs),:)) > 0 
            X   = pSN*data(meg_ch(goodchs),:);
            Xin = X(1:size(SNin,2),:);
            data_in = real(SNin*X(1:size(SNin,2),:));
            if Lout>0
                data_out = real(SNout*X(size(SNin,2)+1:end,:));
            else
                data_out = 0;
            end
            
            data_res = data(meg_ch,:) - (data_in+data_out);
            
            if tsss
                %
                % Temporal analysis
                %
                Bin  = data_in/norm(data_in);
                Bres = data_res/norm(data_res);
                Ein = orth(Bin');
                Eres = orth(Bres');
                [QA,ignore] = qr(Ein,0);
                [QB,ignore] = qr(Eres,0);
                [U,SS,V] = svd(QA'*QB);
                %Up = QA*U';
                Vp = QB*V;
                s = diag(SS);
                inter_indices = find(s>corr_limit);
                length(inter_indices)
                Eproj = Vp(:,inter_indices);
                P = Eproj*Eproj';
                %Bp = ((eye(size(P,1))-P)*data_in')';
                Xp = ((eye(size(P,1))-P)*Xin')';
            else
                %Bp = data_in;
                Xp = Xin;
            end
            
            for k = 1:size(SNin,2)
                Xp(k,:) = Xp(k,:)/norm(Sin(:,k));
            end
                                    
            out_data = real(Sin*Xp);
            out_data(out_mags,:) = out_data(out_mags,:)/(mag2SI*magscale);
            out_data(out_grads,:) = out_data(out_grads,:)/grad2SI;
        else
            out_data = zeros(length(out_meg_ch), size(data, 2));
        end
    end
    
    if ismember(j, Ibar), spm_progress_bar('Set', j); end
    
    spm_progress_bar('Set','ylabel','writing...');
    
    if cont
        Dout(out_meg_ch, from:to)   = out_data;
        Dout(out_other_ch, from:to) = data(other_ch, :);
    else
        Dout(out_meg_ch,  :, j)   = out_data;
        Dout(out_other_ch,  :, j) = data(other_ch, :);
    end
end

%%
Dout = chanlabels(Dout, out_meg_ch, chanlabels(Dref, ref.meg_ch));
Dout = chantype(Dout, out_meg_ch, chantype(Dref, ref.meg_ch));
Dout = units(Dout, out_meg_ch, units(Dref, ref.meg_ch));
Dout = coor2D(Dout, out_meg_ch,  coor2D(Dref, ref.meg_ch));

Dout = chanlabels(Dout, out_other_ch, chanlabels(D, other_ch));
Dout = chantype(Dout,  out_other_ch, chantype(D, other_ch));
Dout = units(Dout, out_other_ch, units(D, other_ch));
Dout = badchannels(Dout,  out_other_ch, badchannels(D, other_ch));
Dout = coor2D(Dout, out_other_ch, coor2D(D, other_ch));

if realign
    SSS = ref.SSS;
    Dout =   sensors(Dout, 'MEG', ref.grad);
    Dout =   fiducials(Dout, ref.fid);
end

[SN_new, sss_indices , nmodes] = basis_condition_adjustment(SSS(1).SN, size(SSS(1).SNin, 2), S.cond_threshold);
pSN = pinv(SN_new);
pSN = pSN(1:nmodes, :);

labelnew = {};
for i = 1:nmodes
    labelnew{i} = ['moment' num2str(i)];
end

mont = [];
mont.labelorg = Dout.chanlabels(out_meg_ch);
mont.labelnew = labelnew;
mont.tra = pSN;

Dout = montage(Dout, 'add', mont);
Dout = chantype(Dout, ':', 'MEG');
Dout = units(Dout, ':', 'fT');

Dout.SSS = SSS;

if ~S.xspace
   Dout = montage(Dout, 'switch', 0); 
end

%-Save the M/EEG dataset
%--------------------------------------------------------------------------
Dout = Dout.history(mfilename, S);

save(Dout);

D = Dout;
%%
%-Cleanup
%--------------------------------------------------------------------------
spm('FigName','TSSS tool: done'); spm('Pointer','Arrow');
spm_progress_bar('Clear');





