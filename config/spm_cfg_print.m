function cfg_print = spm_cfg_print
% SPM Configuration file for 'Print figure'
%__________________________________________________________________________
% Copyright (C) 2005-2012 Wellcome Trust Centre for Neuroimaging

% $Id: spm_cfg_print.m 6924 2016-11-09 11:38:00Z guillaume $

%--------------------------------------------------------------------------
% opts Printing Format
%--------------------------------------------------------------------------
opts         = cfg_menu;
opts.tag     = 'opts';
opts.name    = 'Printing Format';
opts.help    = {['Select the printing format you want. PostScript (PS) is '...
                'the only format that allows to append figures to the same ' ...
                'file.']};
pf           = spm_print('format');
for i=1:numel(pf)
    opts.labels{i} = pf(i).name;
    opts.values{i} = pf(i).label{1};
end
opts.def = @(val)spm_get_defaults('ui.print', val{:});

%--------------------------------------------------------------------------
% fname Print Filename
%--------------------------------------------------------------------------
fname         = cfg_entry;
fname.tag     = 'fname';
fname.name    = 'Print Filename';
fname.strtype = 's';
fname.val     = {''};
fname.help    = {['Filename to print to. If set to an empty string, the ' ...
                'figure will be printed to a file named spm_*.*, in the ' ...
                'current directory. PostScript files will be appended to, ' ...
                'but other files will have ''page numbers'' appended to them.']};

%--------------------------------------------------------------------------
% figname Figure Name
%--------------------------------------------------------------------------
figname           = cfg_entry;
figname.tag       = 'figname';
figname.name      = 'Figure Name';
figname.strtype   = 's';
figname.val       = {'Graphics'};
figname.help      = {['Figure to print. The value entered here must be the ' ...
                    '''Tag'' property of a figure.']};

%--------------------------------------------------------------------------
% fighandle Figure Handle
%--------------------------------------------------------------------------
fighandle         = cfg_entry;
fighandle.tag     = 'fighandle';
fighandle.name    = 'Figure Handle';
fighandle.strtype = 'r';
fighandle.val     = {NaN};
fighandle.help    = {['Figure to print. The value entered here must be a ' ...
                    'figure handle. If it is a non-finite value (Inf/NaN), ' ...
                    'the SPM Graphics window is printed.']};

%--------------------------------------------------------------------------
% fig Figure to print
%--------------------------------------------------------------------------
fig               = cfg_choice;
fig.tag           = 'fig';
fig.name          = 'Figure to print';
fig.values        = {figname, fighandle};
fig.val           = {figname};
fig.help          = {'Figure to print'};

%--------------------------------------------------------------------------
% cfg_print Print figure
%--------------------------------------------------------------------------
cfg_print         = cfg_exbranch;
cfg_print.tag     = 'print';
cfg_print.name    = 'Print figure';
cfg_print.val     = {fname, fig, opts};
cfg_print.prog    = @spm_print;
cfg_print.help    = {'Print figure.'};
