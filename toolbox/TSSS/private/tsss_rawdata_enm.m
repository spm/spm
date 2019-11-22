% tsss_rawdata_enm(infile,outfile,r_sphere,t_window,corr_limit,Lin,Lout,badchs)
% 
% Perform tSSS on rawdata file 'infile' acquired with an Elekta Neuromag
% 306-channel MEG system. The tSSS-processed data are written into 'outfile'.
% The SSS operation is performed in the head coordinate system with the
% expansion origin given by the 3x1 dimensional vector 'r_sphere' ([x y z]
% m; typically [0 0 0.04]). The temporal correlation analysis of tSSS is based 
% on raw data segments of length 't_window' (in seconds) and correlation limit
% 'corr_limit'. The order values of the internal and external SSS bases are 
% 'Lin' and 'Lout', typically 8 and 3, respectively. The list of bad
% channels can be given in the variable 'badchs' in the form 
% badchs = {'MEG XXXX' 'MEG YYYY' 'MEG ZZZZ' ...}.
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
function tsss_rawdata_enm(infile,outfile,r_sphere,t_window,corr_limit,Lin,Lout,badchs)

mag_types = [3022 3023 3024];
grad_types = [3012 3013];
nsamp_max_cov = 1e5;  
magscale = 100;

raw = fiff_setup_read_raw(infile,1);
nsamp = raw.last_samp - raw.first_samp;
data_length = nsamp/raw.info.sfreq;
nwin = round(data_length/t_window);
nsamp_win = nsamp/nwin;  % Number of samples in one window
%
% Extract the channel information
%
if ~isempty(badchs)
    nbad = length(badchs);
else
    nbad = 0;
end
meg_chs = []; mags = []; grads = []; goodchs = [];
for j = 1:raw.info.nchan     
    if raw.info.chs(j).kind == 1;
        goodch = 1;
        meg_chs = [meg_chs j];
        for k = 1:nbad
            if strcmp(raw.info.chs(j).ch_name,badchs{k})
                goodch = 0;
                disp(sprintf('Ignoring channel %s\n',raw.info.chs(j).ch_name));
            end
        end
        if goodch
            goodchs = [goodchs length(meg_chs)];
        end
        if find(raw.info.chs(j).coil_type==grad_types) 
            grads = [grads length(meg_chs)];
        elseif find(raw.info.chs(j).coil_type==mag_types)
            mags = [mags length(meg_chs)];
        end
    end
end
coils = zeros(length(meg_chs));
coils(mags) = 1;
%
% Calculate SSS basis
%
[R,EX,EY,EZ] = fiff_getpos(infile,'head');
disp('Calculating the SSS basis...');
Sin = Sin_vsh_fast(r_sphere,R,EX,EY,EZ,coils,Lin);
Sout = Sout_vsh_fast(r_sphere,R,EX,EY,EZ,coils,Lout);
Sin(mags,:) = magscale*Sin(mags,:);
Sout(mags,:) = magscale*Sout(mags,:);
for j = 1:size(Sin,2) 
    SNin(:,j) = Sin(:,j)/norm(Sin(:,j));
end
for j = 1:size(Sout,2) 
    SNout(:,j) = Sout(:,j)/norm(Sout(:,j));
end
SN = [SNin SNout];
%cond(SN)
pSN = pinv(SN(goodchs,:));
%
% Process the raw data
%
[outfid,cals] = fiff_start_writing_raw(outfile,raw.info);
for j = 1:nwin
    from = (j-1)*nsamp_win + raw.first_samp;
    to = j*nsamp_win - 1 + raw.first_samp;
    data = fiff_read_raw_segment(raw,from,to);
    data(meg_chs(mags),:) = magscale*data(meg_chs(mags),:);
    if size(data,2) >= nsamp_win
        if norm(data) > 0    
            X = pSN*data(meg_chs(goodchs),:);
            data_in = real(SNin*X(1:size(SNin,2),:));
            data_out = real(SNout*X(size(SNin,2)+1:end,:));
            data_res = data(meg_chs,:) - (data_in+data_out);
            %
            % Temporal analysis
            %
            Bin = data_in/norm(data_in);
            Bres = data_res/norm(data_res);
            Ein = orth(Bin');
            Eres = orth(Bres');
            [QA,ignore] = qr(Ein,0);
            [QB,ignore] = qr(Eres,0);
            [U,S,V] = svd(QA'*QB);
            Up = QA*U';
            Vp = QB*V;
            s = diag(S);
            inter_indices = find(s>corr_limit); 
            length(inter_indices)
            Eproj = Vp(:,inter_indices);
            P = Eproj*Eproj';
            Bp = ((eye(size(P,1))-P)*data_in')';
            data(meg_chs,:) = Bp;
            data(mags,:) = (1/magscale)*data(mags,:);
        end
    end
    fiff_write_raw_buffer(outfid,data,cals);
end  
fiff_finish_writing_raw(outfid);







