function con = spm_config_contrasts
% Configuration file for contrast jobs
%_______________________________________________________________________
% Darren Gitelman, baseline 05/02/06
% JA recommended removing extra choice node. Adjusted setupcon code
% accordingly. 05/02/07.
% DRG CS-RCS: $Id: spm_config_contrasts.m,v 1.1 2005-02-08 21:07:14-06 drg Exp drg $

%_______________________________________________________________________

w = spm_jobman('HelpWidth');

%------------------------------------------------------------------------

spm.type = 'files';
spm.name = 'Select SPM.mat';
spm.tag  = 'spmmat';
spm.num  = [1 1];
spm.filter = 'mat';
spm.help   = spm_justify(w,...
    'Select SPM.mat file for contrasts');

name.type    = 'entry';
name.name    = 'Name';
name.tag     = 'name';
name.strtype = 's';
name.num     = [1 1];
name.help    = {'Name of contrast'};

tconvec.type    = 'entry';
tconvec.name    = 'T contrast vector';
tconvec.tag     = 'tconvec';
tconvec.strtype = 's';
tconvec.num     = [1 1];
tconvec.help    = spm_justify(w,...
    'Enter T contrast vector. This is done similarly to the',...
    'SPM2 contrast manager. A 1 x n vector should be entered',...
    'for T-contrasts.');

fconvec.type    = 'entry';
fconvec.name    = 'F contrast vector';
fconvec.tag     = 'fconvec';
fconvec.strtype = 's+';
fconvec.num     = [1 1];
fconvec.help    = spm_justify(w,...
    'Enter F contrast vector. This is done similarly to the',...
    'SPM2 contrast manager. One or multiline contrasts',...
    'may be entered.');

tcon.type   = 'branch';
tcon.name   = 'T-contrast';
tcon.tag    = 'tcon';
tcon.val    = {name,tconvec,};
tcon.help  = {'Enter a T contrast.'};

fcon.type   = 'branch';
fcon.name   = 'F-contrast';
fcon.tag    = 'fcon';
fcon.val    = {name,fconvec,};
fcon.help  = {'Enter a T contrast.'};

consess.type = 'repeat';
consess.name = 'Contrast Sessions';
consess.tag  = 'consess';
consess.values  = {tcon,fcon};
consess.help = {'contrast'};

con.type = 'branch';
con.name = 'Contrast Manager';
con.tag  = 'con';
con.val = {spm,consess};
con.prog   = @setupcon;
con.help = {...
    'Set up T and F contrasts.'};

%--------------------------------------------------------
function setupcon(varargin)

job = varargin{1};

% Change to the analysis directory
%------------------------------------------------------------
if ~isempty(job)
    try
        [pth fn] = fileparts(job.spmmat{:});
    cd(char(pth));
    fprintf('   Changing directory to: %s\n',char(pth));
    catch
        error('Failed to change directory. Aborting contrast setup.')
    end
end

% Load SPM.mat file
%------------------------------------------------------------
load(job.spmmat{:});

for i = 1:length(job.consess)
    if isfield(job.consess{i},'tcon')
    name = job.consess{i}.tcon.name;
    STAT = 'T';
    con  = job.consess{i}.tcon.tconvec(:)';
    else %fcon
        name = job.consess{i}.fcon.name;
        STAT = 'F';
        con = cellstr(job.consess{i}.fcon.fconvec);
        
    end
    
    % Basic checking of contrast
    %------------------------------------------------------------
    [c,I,emsg,imsg] = spm_conman('ParseCon',con,SPM.xX.xKXs,STAT);
    
    % Fill-in the contrast structure
    %------------------------------------------------------------
    if all(I)
		DxCon = spm_FcUtil('Set',name,STAT,'c',c,SPM.xX.xKXs);
	else
		DxCon = [];
    end
    
    % Append to SPM.xCon. SPM will automatically save any contrasts that
    % evaluate successfully.
    %------------------------------------------------------------
    if isempty(SPM.xCon),
        SPM.xCon = DxCon;
    else
        SPM.xCon(end+1) = DxCon;
    end;
    spm_contrasts(SPM,length(SPM.xCon));
    
end


