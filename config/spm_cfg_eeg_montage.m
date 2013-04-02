function S = spm_cfg_eeg_montage
% configuration file for reading montage files
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_cfg_eeg_montage.m 5377 2013-04-02 17:07:57Z vladimir $

D = cfg_files;
D.tag = 'D';
D.name = 'File Name';
D.filter = 'mat';
D.num = [1 1];
D.help = {'Select the EEG mat file.'};

montagefile = cfg_files;
montagefile.tag = 'montagefile';
montagefile.name = 'Montage file name';
montagefile.filter = 'mat';
montagefile.num = [1 1];
montagefile.help = {'Select a montage file.'};

keepothers = cfg_menu;
keepothers.tag = 'keepothers';
keepothers.name = 'Keep other channels';
keepothers.labels = {'Yes', 'No'};
keepothers.values = {1, 0};
keepothers.val = {0};
keepothers.help = {'Specify whether you want to keep channels that are not contributing to the new channels'};

montage = cfg_branch;
montage.tag = 'montage';
montage.name = 'Load montage from file';
montage.val = {montagefile, keepothers};

montagename         = cfg_entry;
montagename.tag     = 'montagename';
montagename.name    = 'Montage name';
montagename.help    = {'Specify the name of existing montage'};
montagename.strtype = 's';
montagename.num     = [1 Inf];

montagenum         = cfg_entry;
montagenum.tag     = 'montagenum';
montagenum.name    = 'Montage index';
montagenum.help    = {'Specify the index of existing montage'};
montagenum.strtype = 'n';
montagenum.num     = [1 1];

montspec = cfg_choice;
montspec.tag = 'montspec';
montspec.name = 'Montage specification';
montspec.values = {montage, montagename, montagenum};
montspec.val = {montage};
montspec.help = {'Choose how to specify the montage'};

blocksize = cfg_entry;
blocksize.tag = 'blocksize';
blocksize.name = 'Block size';
blocksize.strtype = 'r';
blocksize.val = {655360};
blocksize.num = [1 1];
blocksize.help = {'size of blocks used internally to split large files default ~20Mb'};

prefix         = cfg_entry;
prefix.tag     = 'prefix';
prefix.name    = 'Filename Prefix';
prefix.help    = {'Specify the string to be prepended to the filenames of the montaged ',...
    'dataset (if generated). Default prefix is ''M''.'};
prefix.strtype = 's';
prefix.num     = [1 Inf];
prefix.val     = {'M'};

write = cfg_branch;
write.tag = 'write';
write.name = 'Write';
write.val = {montspec, blocksize, prefix};
write.help = {'Apply montage and write out to a new dataset.'};

online = cfg_branch;
online.tag = 'switch';
online.name = 'Switch';
online.val = {montspec};
online.help = {'Switch to new online montage without writing out a new file'};

add = cfg_branch;
add.tag = 'add';
add.name = 'Add';
add.val = {montage};
add.help = {'Add montage to the set of online montages without applying'};

clr = cfg_const;
clr.tag = 'clear';
clr.name = 'Clear';
clr.val  = {1};
clr.help = {'Clear all online montages and revert to the original state.'};

mode = cfg_choice;
mode.tag = 'mode';
mode.name = 'Mode';
mode.values = {write, online, add, clr};
mode.val = {write};
mode.help = {'Choose between writing a new dataset or online montage operation'};

S = cfg_exbranch;
S.tag = 'montage';
S.name = 'Montage';
S.val = {D, mode};
S.help = {'Apply a montage (linear transformation) to EEG/MEG data.'};
S.prog = @eeg_montage;
S.vout = @vout_eeg_montage;
S.modality = {'EEG'};

function out = eeg_montage(job)
% construct the S struct
S.D = job.D{1};
S.mode = char(fieldnames(job.mode));
switch S.mode
    case 'write'
        S.blocksize = job.mode.write.blocksize;
        S.prefix = job.mode.write.prefix;
        switch char(fieldnames(job.mode.write.montspec))
            case 'montagename'
                S.montage = job.mode.write.montspec.montagename;
            case 'montagenum'
                S.montage = job.mode.write.montspec.montagenum;
            case 'montage'
                S.montage = char(job.mode.write.montspec.montage.montagefile);
                S.keepothers = job.mode.write.montspec.montage.keepothers;
        end
    case 'switch'
        switch char(fieldnames(job.mode.switch.montspec))
            case 'montagename'
                S.montage = job.mode.switch.montspec.montagename;
            case 'montagenum'
                S.montage = job.mode.switch.montspec.montagenum;
            case 'montage'
                S.montage = char(job.mode.switch.montspec.montage.montagefile);
                S.keepothers = job.mode.switch.montspec.montage.keepothers;
        end
    case 'add'
        S.montage = char(job.mode.add.montage.montagefile);
        S.keepothers = job.mode.add.montage.keepothers;
end
out.D = spm_eeg_montage(S);
out.Dfname = {fullfile(out.D.path, out.D.fname)};

function dep = vout_eeg_montage(job)
% Output is always in field "D", no matter how job is structured
dep = cfg_dep;
dep.sname = 'Montaged data';
% reference field "D" from output
dep.src_output = substruct('.','D');
% this can be entered into any evaluated input
dep.tgt_spec   = cfg_findspec({{'strtype','e'}});

dep(2) = cfg_dep;
dep(2).sname = 'Montaged Datafile';
% reference field "Dfname" from output
dep(2).src_output = substruct('.','Dfname');
% this can be entered into any file selector
dep(2).tgt_spec   = cfg_findspec({{'filter','mat'}});


