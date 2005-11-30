function [sdip,fit_opt,Psave] = spm_eeg_inv_ecd_fitDip_ui(D)

% FORMAT [sdip,set_loc_o,Psave] = spm_eeg_inv_ecd_fitDip_ui(D)
%
% 'spm_eeg_inv_ecd_fitDip_ui' provides the GUI for the dipole fitting 
% routine 'spm_eeg_inv_ecd_fitDipS' and 'spm_eeg_inv_ecd_fitDip'. 
% The former uses a "Realistic sphere forward solution", the later 
% a BEM solution (which is not fully supported yet!).
%
% If the data structure is passed then it uses the information it contains,
% then it asks the user to specify (if it's not defined yet):
%   - the EEG data, and the time window to be used (VERY important)
%   - the head model, comprising surfaces, electrodes & conductivities info.
%   - the brain mask to limit the dipole locations
%   - nr of dipoles to be used
%   - nr of random starting locations to be considered
%   - orientation of the sources over the time window (free or fixed)
%
% Output :
%   - sdip: fitted dipole(s) strucure
%       + n_seeds: # seeds set used
%       + n_dip: # fitted dipoles on the EEG time series
%       + loc: location of fitted dipoles, cell{1,n_seeds}(3 x n_dip)
%       + L: leadfiled associated with each set of n_dip fitted dipoles, 
%            cell(1,n_seeds){N_el x 3*n_dip}
%       + j: sources amplitude over the time window, cell{1,n_seeds}(3*n_dip x Ntimebins)
%       + res: residuals of each fit, vect(1,n_seeds)
%       + cost: cost of each fit, vect(1,n_seeds)
%       + rres: relative residuals of each fit, vect(1,n_seeds)
%       + varexpl: variance explained by fitted dipole(s)
%       + M: 4x4 affine transformation matrix mapping from brain image voxel coordinates
%            to real world coordinates
%       + Mtb: index of maximum power in EEG time series used
%       + fit_opt: fitting options used
%       + wind_be: time window defining timeseries used, in complete data set
%       + Pdata: data file
%       + fname: name of the file where the structure is saved.
%   - set_loc_o: set of starting location priors, used as seeds.
%   - Psave: file name of saved file.
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips,
% $Id$

% Loading various bits
%_____________________
spm('Clear')
spm('FigName','Fitting ECD on EEG data') ;

q_model=2; % At the moment, only consider the spherical model
if nargin==0
    % Select the head model
    % [q_model,pos] = spm_input('Model type :','+1','realistic BEM|3 fitted sph ',[1,2],2) ;
    % At the moment, use only the realistic sphere!
    if q_model==1
        Pmod = spm_select(1,'^model.*\ifs.*\.mat$','Model mat file');
    else
        Pmod = spm_select(1,'^model.*\Rs.*\.mat$','Model mat file');
    end
else
    Pmod = D.inv{end}.model;
end
load(Pmod)

if nargin==0
    % Select the data
    [q_data,pos] = spm_input('Data format :','+1','in-house bin|simple mat file',[1,2],2) ;
    if q_data==1
        D = spm_eeg_ldata
    end
else
    q_data = 1
end

if q_data==1
    Pdata = D.fname;
    ss = size(D.data) ;
    DNchan = ss(1);
    Dtb = ss(2);
    if length(ss)==3
        if ss(3)>1
            [cond,pos] = spm_input(['Condition to use: ',num2str(1:ss(3)),' ?'],'+1','e',1) ;
        else
            cond = 1;
        end
    else
        cond = 1;
    end
elseif q_data==2
    Pdata = spm_select(1,'mat','Data mat file');
    load(Pdata);
    [DNchan,Dtb] = size(data);
end

if q_data==1
    % Check the data, channels to use, nr of electrodes, etc
    Use_chan = 1:D.Nchannels;
    chan_to_rem = [];
    if isfield(D.channels,'heog')
        chan_to_rem = [chan_to_rem D.channels.heog];
    end
    if isfield(D.channels,'veog')
        chan_to_rem = [chan_to_rem D.channels.veog];
    end
    if isfield(D.channels,'eog')
        chan_to_rem = [chan_to_rem D.channels.eog];
    end
    if isfield(D.channels,'ecg')
        chan_to_rem = [chan_to_rem D.channels.ecg];
    end
    if isfield(D.channels,'others')
        chan_to_rem = [chan_to_rem D.channels.others];
    end
    if isfield(D.channels,'Bad')
        chan_to_rem = [chan_to_rem D.channels.Bad'];
    end
    if isfield(D.channels,'reference')
        chan_to_rem = [chan_to_rem D.channels.reference];
    end
    Use_chan(chan_to_rem) = []; % and possibly some other channels.
    
    % NUse_chan = length(D.cov.Use_chan);
    NUse_chan = length(Use_chan);
    
    order_dat2mod = zeros(1,model.electrodes.nr); missing_chan = [];
    % Map name of channels in data onto channel order specified in model
    for i = 1:model.electrodes.nr
        index = [];
        for j = 1:D.Nchannels
            if ~isempty(find(strcmpi(D.channels.name{j}, model.electrodes.names{i})))
                index = [index j];
            end
        end
        
        if isempty(index)
            warning(sprintf('No channel named %s found in data.', model.electrodes.names{i}));
            missing_chan = [missing_chan i];
        else
            % take only the first found channel descriptor
            order_dat2mod(i) = index(1);
        end
    end
    %     D.channels.name(order_dat2mod)
elseif q_data==2
    disp('I assume the data channels are in same order as model !')
end

% Use a reduce set of electrodes ?
% dNchannels = model.electrodes.nr-DNchan ;
dNchannels = model.electrodes.nr-NUse_chan ;
if dNchannels<0
    spm('alert*',{'Your data have more channels than the model,' ...
            'I can''t deal with this...'},'Nchannels error')
    return
elseif dNchannels>0
    spm('alert!',{['Your data have less channels (',num2str(NUse_chan),...
                ') than the model (',num2str(model.electrodes.nr),'),'] ...
            [num2str(dNchannels),...
                ' in the leadfield matrix will be ''removed''.']},'Nchannels warning')
    try
        rem_electr = [];
        for ii=chan_to_rem
            p_rem = find(order_dat2mod==ii);
            if ~isempty(p_rem)
                rem_electr = [rem_electr p_rem];
            end
        end
        if length(rem_electr)~=dNchannels
            error('Input channels to remove manually!')
        end
        keep_electr = 1:model.electrodes.nr;
        keep_electr(rem_electr) = [];
    catch
        flag = 0; pos = pos+1;
        while flag==0
            [rem_electr,pos] = spm_input(['List (',num2str(dNchannels), ...
                    ') of electrodes to remove in the model'],pos,'e',missing_chan,dNchannels);
            keep_electr = 1:model.electrodes.nr;
            keep_electr(rem_electr) = [];
            if length(keep_electr)==DNchan, flag=1; end
        end
    end
    order_dat2mod(rem_electr) = [];
    % 	% Modifying the electrodes structure
    model.electrodes.vert(rem_electr) = [];
    model.electrodes.tri(rem_electr) = [];
    model.electrodes.XYZmm(:,rem_electr) = [];
    
    model.electrodes.names(rem_electr,:) = [];
    model.electrodes.nr = NUse_chan;
    model.electrodes.info = [model.electrodes.info,' ; some electrodes removed ''cos of data'];
    
    if q_model==1
        % Modifying the IFS matrices
        tr_rem_electr = [3*rem_electr-2 ; 3*rem_electr-1 ; 3*rem_electr];
        tr_rem_electr = tr_rem_electr(:);
        for ii=1:length(model.IFS)
            model.IFS{ii}(tr_rem_electr,:) = [];
        end
    else
        model.spheres.Sc_elXYZ(:,rem_electr) = [];
    end
    if q_data==1
        order_dat2mod = order_dat2mod(find(order_dat2mod));
    end
else
    rem_electr = [];
end

% Window of time used
flag=0; pos = pos+1;
while flag==0
    if q_data==1
        % Express time window in ms.
        ms_tb = 1000/D.Radc ; % milisecond per time bin
        b_ms = -D.events.start*ms_tb;
        e_ms = D.events.stop*ms_tb;
        [wind_be_ms,pos] = spm_input(...
                ['Time window to use (in ms, [',num2str(ms_tb),'])'], ...
                pos,'e',round([b_ms e_ms]));
        wind_be = round((wind_be_ms-b_ms)/ms_tb+1);
    else
        [wind_be,pos] = spm_input(['Time window to use'],pos,'e',[1 Dtb]);
    end
    if length(wind_be)==1
        if wind_be(1)>0 && wind_be(1)<=Dtb
            wind_be(2) = wind_be(1); 
            flag = 1;
        end
    elseif length(wind_be)==2
        if wind_be(2)>wind_be(1) && all(wind_be>0) && all(wind_be<=Dtb)
            flag = 1;
        end
    end
end

% Number of dipoles per set
flag=0; pos = pos+1;
while flag==0
    [n_dip,pos] = spm_input(['Number of dipoles to be used'],pos,'e',1);
    if n_dip<=floor(DNchan/6), flag=1; end
    if (n_dip<=floor(DNchan/6)) && (n_dip>floor(DNchan/12)),
        spm('alert!',{['You really have a lot of dipoles to fit (',num2str(n_dip),')'] ...
            ['compared to the number of channels (',num2str(DNchan),')'] ...
            ['Results may be unstable and unreliable']},'Ndipoles warning')
    end
end

% Number of random seed sets
[n_seeds,pos] = spm_input(['Number of random seeding locations to be used (1 for a priori loc.)'],'+1','e',n_dip*20);

% If only 1 seed, I should define myself the a priori location, and possibly fix these locations.
% In this last case, I should be able to fix the orientation too.
q_fxd_loc = 0; q_fxd_or = 0; set_loc_o = []; fxd_or = [];
if n_seeds==1
    set_loc_o = zeros(3,n_dip);
    %     set_loc_1seed = [[0 0 0]' [-30 0 0]' [30 0 0]' [0 0 30]' [0 30 0]' [0 -30 0]' ];
    set_loc_1seed = [[-.3 .9 42.8]' [41 -4 47]' [30 0 0]' [0 0 30]' [0 30 0]' [0 -30 0]' ];
    for ii=1:n_dip
        [set_loc_o(:,ii)] = spm_input(['A priori location of dipole # ',num2str(ii)],pos+1, ...
            'e',set_loc_1seed(:,ii)',3);
    end
    pos = pos+1;
    [q_fxd_loc,pos] = spm_input('Leave the dipoles location fixed ?','+1','y/n',[1,0],2);
    if q_fxd_loc
        [q_fxd_or,pos] = spm_input('Fix the dipoles orientation ?','+1','y/n',[1,0],2);
        if q_fxd_or % Define the dipoles orientation
            fxd_or = zeros(3,n_dip);
            for ii=1:n_dip
                tmp_or = spm_input(['A priori orientation of dipole # ',num2str(ii)],pos+1,'e',[],3);
                fxd_or(:,ii) = tmp_or/norm(tmp_or);
            end
            pos = pos+1;
        end
    end
end

% Orientation of the dipoles
if (wind_be(2)-wind_be(1))>1 && ~q_fxd_or
    text_opt = ['totally free orientation|',...
            'fixed orientation, weighted (EEG pow) mean of orientations over time|',...
            'fixed orientation, as at the max of EEG power'];
    [or_opt,pos] = spm_input('Dipoles orientation over time?',pos+1,'m',text_opt);
elseif q_fxd_or
    or_opt = 4;
else
    or_opt = 1;
end

% Select obrain volume to constrain sources
if nargin==0
    Vbr = spm_vol(spm_select(1,'^.*obrain.*\..*$','Brain mask image'));
else
    Vbr = spm_vol(D.inv{end}.mesh.msk_iskull);
    if isempty(Vbr)
        Vbr = spm_vol(model.flags.fl_bin.Pisk);
    end
end

Psdip = ['S',num2str(n_dip),'dip_',spm_str_manip(Pdata,'rt'), ...
        '_n',num2str(n_seeds),'_o',num2str(or_opt),'_t',num2str(wind_be(1)),'_',num2str(wind_be(2))];
[Psave,pos] = spm_input(['Save results in'],'+1','s',Psdip);

if ~isfield(model.param,'sigma')
    model.param.sigma = spm_input('Conductivity for the br-sk-sc volumes','+1','e',[.33 .004 .33])
end

% Fitting the dipoles !
%______________________
% A few options passed into one simple structure.
fit_opt.q_model = q_model;          % Model used: BEM (1), analytical 3 fitted sph (2)
fit_opt.or_opt = or_opt;            % Def: 1, free orientation
fit_opt.rem_electr = rem_electr;    % Def: none
fit_opt.q_fxd_loc = q_fxd_loc;      % Def: 0
fit_opt.set_loc_o = set_loc_o;      % Def: none
fit_opt.q_fxd_or = q_fxd_or;        % Def: 0
fit_opt.fxd_or = fxd_or;            % Def: none

if q_model==1
    if q_data==1
        [sdip,fit_opt] = ...
            spm_eeg_inv_ecd_fitDip(squeeze(D.data(order_dat2mod, ...
                wind_be(1):wind_be(2),cond(1))), ...
                model,Vbr,n_dip,n_seeds,fit_opt);
    elseif q_data==2
        [sdip,fit_opt] = spm_eeg_inv_ecd_fitDip(data(:,wind_be(1):wind_be(2)), ...
            model,Vbr,n_dip,n_seeds,fit_opt);
    end  
elseif q_model==2
    if q_data==1
        [sdip,fit_opt] = ...
            spm_eeg_inv_ecd_fitDipS(squeeze(D.data(order_dat2mod, ...
                                wind_be(1):wind_be(2),cond(1))), ...
                                model,Vbr,n_dip,n_seeds,fit_opt);
    elseif q_data==2
        [sdip,fit_opt] = spm_eeg_inv_ecd_fitDipS(data(:,wind_be(1):wind_be(2)), ...
            model,Vbr,n_dip,n_seeds,fit_opt);
    end
end
sdip.wind_be = wind_be;
sdip.Pdata = Pdata;
sdip.fname = Psave;

% Save results
%______________
save(Psave,'sdip','fit_opt')

