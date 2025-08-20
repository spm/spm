function bsd_cfg = spm_bsd_batch
% SPM Configuration file for Bayesian Spectral Decomposition (BSD)
%__________________________________________________________________________

% Directory and output name
%--------------------------------------------------------------------------
% dir Directory
datafiles         = cfg_files; 
datafiles.tag     = 'Dfiles'; 
datafiles.name    = 'Dataset(s)';
datafiles.help    = {'Select MEEG dataset files.'};
datafiles.filter  = 'mat';
datafiles.ufilter = '.*';
datafiles.num     = [1 Inf];

modality        = cfg_menu; 
modality.tag    = 'modality'; 
modality.name   = 'Modality'; 
modality.labels = {'EEG', 'MEG', 'MEGPLANAR', 'LFP'}; 
modality.values = {'EEG', 'MEG', 'MEGPLANAR', 'LFP'}; 

% Tdcm Time window
Tdcm         = cfg_entry;
Tdcm.tag     = 'Tdcm';
Tdcm.name    = 'Time window (ms)';
Tdcm.help    = {'Specify the time window in milliseconds as [start end].'};
Tdcm.strtype = 'r';
Tdcm.num     = [0 Inf];
Tdcm.val     = {[0 0]}; 

% Fdcm Frequency window
Fdcm         = cfg_entry;
Fdcm.tag     = 'Fdcm';
Fdcm.name    = 'Frequencies (Hz)';
Fdcm.help    = {'Specify the frequencies in Hz.'};
Fdcm.strtype = 'r';
Fdcm.num     = [0 Inf];
Fdcm.val     = {1:64}; 

chind           = cfg_entry;
chind.tag       = 'chind'; 
chind.name      = 'Specific channel indices'; 
chind.val       = {[]};
chind.strtype   = 'n';
chind.num       = [0 Inf]; 

chnames      = cfg_entry; 
chnames.tag  = 'chnames'; 
chnames.name = 'Specific channel labels';
chnames.val  = {};
chnames.strtype = 's+'; 
chnames.num  = [0, Inf]; 

chall        = cfg_const; 
chall.tag    = 'chall'; 
chall.name   = 'All channels'; 
chall.val    = {'all'};

chsel        = cfg_choice;
chsel.tag    = 'chsel';
chsel.name   = 'Channel selection'; 
chsel.values = {chall, chnames, chind};
chsel.val    = {chall}; 

% condidx Condition indices
condidx         = cfg_entry;
condidx.tag     = 'trials'; 
condidx.name    = 'Condition indices'; 
condidx.strtype = 'n'; 
condidx.num     = [0 Inf];
condidx.val     = {1};
condidx.help    = {'Index of the conditions to model.'}; 

% D Time bin decimation
D         = cfg_entry;
D.tag     = 'D';
D.name    = 'Time decimation';
D.help    = {'Resample by decimating by a factor (e.g., 1 for no decimation, 2 for half).'};
D.strtype = 'n';
D.val     = {1};
D.num     = [1 1];

datameeg      = cfg_branch; 
datameeg.tag  = 'meeg'; 
datameeg.name = 'MEEG dataset(s)'; 
condidx.help    = {'Index of the conditions to model.'}; 

datameeg.val  = {datafiles, Tdcm, Fdcm, modality, chsel, condidx, D}; 

arrayvals         = cfg_entry; 
arrayvals.tag     = 'y';
arrayvals.name    = 'Power spectra';
arrayvals.help    = {'Specify the data matrix.'};
arrayvals.strtype = 'n';
arrayvals.num     = [Inf Inf];

arrayfreq         = cfg_entry; 
arrayfreq.tag     = 'Hz';
arrayfreq.name    = 'Frequency vector';
arrayfreq.help    = {'Specify the frequency vector.'};
arrayfreq.strtype = 'n';
arrayfreq.num     = [1 Inf];

dataarray      = cfg_branch; 
dataarray.tag  = 'xY';
dataarray.name = 'Data array'; 
dataarray.val  = {arrayvals, arrayfreq}; 

datasource        = cfg_choice; 
datasource.tag    = 'source'; 
datasource.name   = 'Data source'; 
datasource.values = {datameeg, dataarray};
datasource.val    = {datameeg}; 

%--------------------------------------------------------------------------

fqcustom         = cfg_entry; 
fqcustom.tag     = 'val';
fqcustom.name    = 'Frequency range';
fqcustom.strtype = 'n';
fqcustom.num     = [1 2];

fqdelta = cfg_entry; 
fqdelta.tag = 'val'; 
fqdelta.name = 'Delta (below 4Hz)';
fqdelta.val  = {[0 4]}; 

fqtheta = cfg_const; 
fqtheta.tag = 'val'; 
fqtheta.name = 'Theta (from 4Hz to 8Hz)';
fqtheta.val  = {[4 8]}; 

fqalpha = cfg_const; 
fqalpha.tag = 'val'; 
fqalpha.name = 'Alpha (from 8Hz to 12Hz)';
fqalpha.val  = {[8 12]}; 

fqbeta = cfg_const; 
fqbeta.tag = 'val'; 
fqbeta.name = 'Beta (from 12Hz to 30Hz)';
fqbeta.val  = {[12 30]}; 

fqgamma = cfg_const; 
fqgamma.tag = 'val'; 
fqgamma.name = 'Gamma (from 30Hz to 64Hz)';
fqgamma.val  = {[30 64]};

frequencybands        = cfg_repeat; 
frequencybands.tag    = 'freqs';
frequencybands.name   = 'Frequency bands'; 
frequencybands.values = {fqcustom, fqdelta, fqtheta, fqalpha, fqbeta, fqgamma};
frequencybands.num    = [0 Inf]; 
frequencybands.help   = {['Specify frequency bands in which to fit a ' ...
    'unique peak. Frequency bands cannot overlap. ']}; 

% Condition specific effects
%--------------------------------------------------------------------------
% design Design matrix
design         = cfg_entry;
design.tag     = 'X';
design.name    = 'Design matrix';
design.help    = {'Specify the design matrix.'};
design.strtype = 'n';
design.num     = [Inf Inf];
design.val     = {[]}; 


% condyes Yes
condyes      = cfg_branch;
condyes.tag  = 'condyes'; 
condyes.name = 'Yes';
condyes.val  = {design}; 

% condno No
condno      = cfg_const;
condno.tag  = 'condno';
condno.name = 'No';
condno.val  = {false};
condno.help = {''};

% condchoice Model condition specific effects?
condchoice        = cfg_choice; 
condchoice.tag    = 'condchoice'; 
condchoice.name   = 'Model condition specific effects?';
condchoice.values = {condyes, condno}; 
condchoice.val    = {condno}; 

% Spatial model
%--------------------------------------------------------------------------
Lpos             = cfg_entry; 
Lpos.tag         = 'Lpos'; 
Lpos.name        = 'Source coordinates (mm)'; 
Lpos.help        = {'Specify the source coordinates.'};
Lpos.strtype     = 'r';
Lpos.num         = [Inf 3];

% Sname Source Names
Sname         = cfg_entry;
Sname.tag     = 'Sname';
Sname.name    = 'Source names';
Sname.help    = {'Enter the names of sources as a cell array of strings.'};
Sname.strtype = 's+';
Sname.num     = [0 Inf];

% Nmodes Number of spatial modes
Nmodes         = cfg_entry;
Nmodes.tag     = 'Nmodes';
Nmodes.name    = 'Number of spatial modes';
Nmodes.help    = {['Specify the number of spatial modes to use. One BSD ' ...
    'model will be fit per spatial mode.']};
Nmodes.strtype = 'n';
Nmodes.val     = {8};
Nmodes.num     = [1 1];

spatialch        = cfg_branch; 
spatialch.tag    = 'spatialch'; 
spatialch.name   = 'Channel space'; 
spatialch.val    = {Nmodes};

spatialecd       = cfg_const; 
spatialecd.tag   = 'typelfp';
spatialecd.name  = 'Local Field Potential (LFP)'; 
spatialecd.val   = {'LFP'};

spatialecd       = cfg_const; 
spatialecd.tag   = 'typeecd';
spatialecd.name  = 'Equivalent Current Dipole (ECD)'; 
spatialecd.val   = {'ECD'};

spatialimg        = cfg_const; 
spatialimg.tag    = 'typeimg';
spatialimg.name   = 'Imaging (IMG)'; 
spatialimg.val    = {'IMG'};

modelchoice        = cfg_choice;
modelchoice.tag    = 'model'; 
modelchoice.name   = 'Model type';
modelchoice.values = {spatialimg, spatialecd};
modelchoice.val    = {spatialecd}; 

spatialsrc        = cfg_branch; 
spatialsrc.tag    = 'spatialsrc'; 
spatialsrc.name   = 'Source space (use forward model)'; 
spatialsrc.val    = {modelchoice, Lpos, Sname, Nmodes};

spatialchoice        = cfg_choice;
spatialchoice.tag    = 'spatialchoice'; 
spatialchoice.name   = 'Data space';
spatialchoice.values = {spatialch, spatialsrc};
spatialchoice.val    = {spatialch}; 


% Model options
%--------------------------------------------------------------------------

% separatenull Fit null model separately
separatenull         = cfg_menu;
separatenull.tag     = 'separatenull';
separatenull.name    = 'Initial aperiodic fit';
separatenull.help    = {'Specify whether to initially fit the aperiodic component.'};
separatenull.labels  = {'Yes', 'No'};
separatenull.values  = {true, false};
separatenull.val     = {true};

% fitlog Fit log power spectra
fitlog         = cfg_menu;
fitlog.tag     = 'fitlog';
fitlog.name    = 'Log-transform power spectra';
fitlog.help    = {'Specify whether to fit log-transformed power spectra.'};
fitlog.labels  = {'Yes', 'No'};
fitlog.values  = {true, false};
fitlog.val     = {true};

% powerline Filter for power line noise
powerline         = cfg_entry;
powerline.tag     = 'powerline';
powerline.name    = 'Powerline filter (Hz)';
powerline.help    = {'Specify the range for the power line filter as [start end] in Hz.'};
powerline.strtype = 'n';
powerline.val     = {[49 51]}; 
powerline.num     = [1 2];

% fitlog Fit log power spectra
plot         = cfg_menu;
plot.tag     = 'plot';
plot.name    = 'Plot';
plot.help    = {'Specify whether to plot while fitting.'};
plot.labels  = {'Yes', 'No'};
plot.values  = {true, false};
plot.val     = {true};

% fitlog Fit log power spectra
print         = cfg_menu;
print.tag     = 'print';
print.name    = 'Print';
print.help    = {'Specify whether to print while fitting.'};
print.labels  = {'Yes', 'No'};
print.values  = {true, false};
print.val     = {true};


options     = cfg_branch;
options.tag = 'options';
options.name = 'Model options';
options.val = {condchoice, spatialchoice, separatenull, fitlog, powerline, plot, print}; 

% prefix Output prefix
prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Output prefix';
prefix.help    = {'Specify a prefix for the output.'};
prefix.strtype = 's';
prefix.num     = [0 Inf];
prefix.val     = {'BSD_'};


% Create configuration branch
%--------------------------------------------------------------------------
bsd_cfg      = cfg_exbranch;
bsd_cfg.tag  = 'spm_bsd';
bsd_cfg.name = 'BSD - Specify and estimate';
bsd_cfg.val  = {datasource, frequencybands, options, prefix};
bsd_cfg.help = {'Configure and run Bayesian Spectral Decomposition (BSD).'};
bsd_cfg.prog = @spm_run_bsd;
bsd_cfg.vout = @vout_spm_bsd;

end

% Output function
%--------------------------------------------------------------------------
function vout = vout_spm_bsd(job)
vout = cfg_dep;
vout.sname = 'BSD Output';
vout.src_output = substruct('.','BSD');
vout.tgt_spec = cfg_findspec({{'filter','mat','strtype','e'}});
end

function BSD = spm_run_bsd(job)
    disp(job)

    if isfield(job.source, 'meeg') && numel(job.source.meeg.Dfiles) > 1
        BSD = cell(numel(job.source.meeg.Dfiles), 1);
        for i = 1:numel(job.source.meeg.Dfiles)
            job_ = job;
            job_.source.meeg.Dfiles = job.source.meeg.Dfiles(i); 
            BSD{i} = spm_run_bsd(job); 
        end
        return 
    end
    
    BSD = []; 

    if isfield(job.source, 'meeg') 
        BSD.xY.Dfile       = job.source.meeg.Dfiles{1}; 
        BSD.xY.modality    = job.source.meeg.modality; 
        BSD.options.trials =  job.source.meeg.trials; 
        BSD.options.D      =  job.source.meeg.D; 
        BSD.options.Tdcm   =  job.source.meeg.Tdcm; 
        BSD.options.Fdcm   =  job.source.meeg.Fdcm; 


        if isfield(job.source.meeg.chsel, 'chind')
            BSD.xY.Ic = job.source.meeg.chsel.chind; 
        elseif isfield(job.source.meeg.chsel, 'chnames')
            BSD.xY.name = cellstr(job.source.meeg.chsel.chnames); 
        end
 
        BSD.name = spm_file(BSD.xY.Dfile, 'prefix', job.prefix);   
    else 
        BSD.xY.Hz = job.source.xY.Hz;
        BSD.xY.y  = job.source.xY.y;
    end

    BSD.fqs = cellfun(@(c) c.val, job.freqs, 'UniformOutput', false); 

    if isfield(job.options.spatialchoice, 'spatialch')
        BSD.options.spatial = 'Chan'; 
    else
        BSD.Lpos         = job.options.spatialchoice.spatialsrc.Lpos;
        BSD.Sname        = job.options.spatialchoice.spatialsrc.Sname;
        
        if isfield( job.options.spatialchoice.spatialsrc.model, 'typelfp')
            BSD.options.spatial = 'LFP';
        elseif isfield( job.options.spatialchoice.spatialsrc.model, 'typeecd')
            BSD.options.spatial = 'ECD';
        else
            BSD.options.spatial = 'IMG';
        end
        BSD.options.Nmodes = job.options.spatialchoice.spatialsrc.Nmodes;
    end

    if isfield(job.options.condchoice, 'condyes')
        BSD.xU.X           = job.options.condchoice.condyes.X; 
    end

    BSD.options.fitlog       = job.options.fitlog;
    BSD.options.powerline    = job.options.powerline;
    BSD.options.separatenull = job.options.separatenull;
    BSD.options.noprint      = ~job.options.print;
    BSD.options.nograph      = ~job.options.plot;

    BSD = spm_bsd(BSD); 
end