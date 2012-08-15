function vdm = FieldMap_Run(job)
% Auxillary file for running FieldMap jobs
%
% FORMAT vdm = Fieldmap_Run(job)
%
% job  - FieldMap job structure containing following elements
%        common to all jobs and specific to the type of job
%     Common to all jobs:
%        defaults - cell array containing name string of the defaults file
%        options - structure containing the following:
%           epi - cell array containing name string of epi image to unwarp
%           matchvdm - match vdm to epi or not (1/0)
%           writeunwarped - write unwarped EPI or not (1/0)
%           anat - cell array containing name string of anatomical image
%           matchanat - match anatomical image to EPI or not (1/0)
%
%     Elements specific to job type:
%        precalcfieldmap - name of precalculated fieldmap
%
%        phase - name of phase image for presubtracted phase/mag job
%        magnitude - name of magnitude image for presubtracted phase/mag job
%
%        shortphase - name of short phase image for phase/mag pair job
%        longphase - name of short phase image for phase/mag pair job
%        shortmag - name of short magnitude image for phase/mag pair job
%        longmag - name of short magnitude image for phase/mag pair job
%
%        shortreal - name of short real image for real/imaginary job
%        longreal - name of long real image for real/imaginary job
%        shortimag - name of short imaginary image for real/imaginary job
%        longimag - name of long imaginary image for real/imaginary job
%__________________________________________________________________________
% Copyright (C) 2007 Wellcome Trust Centre for Neuroimaging

% Chloe Hutton & Jesper Andersson
% $Id: FieldMap_Run.m 4844 2012-08-15 19:12:38Z guillaume $


%--------------------------------------------------------------------------
% Set up default parameters and structures
%--------------------------------------------------------------------------

% Open the FieldMap control window with visibility off. This allows the
% graphics display to work.
FieldMap('Welcome','Off');
%IP = FieldMap('Initialise'); % Gets default params from pm_defaults

% Here load the selected defaults file if selected
if isfield(job.defaults,'defaultsfile')
    m_file = job.defaults.defaultsfile;
    m_file = spm_file(m_file,'cpath');
    %m_file = spm_file(m_file{1},'ext','');
    %m_file = spm_file(m_file{1},'basename');
    pm_defs = FieldMap('SetParams',m_file); % Gets default params from pm_defaults
elseif isfield(job.defaults,'defaultsval')
    pm_defs = job.defaults.defaultsval;
    echotimes=pm_defs.et;
    pm_defs.et=[];
    pm_defs.et{1}=echotimes(1);
    pm_defs.et{2}=echotimes(2);
    pm_defs.uflags.etd=pm_defs.et{2}-pm_defs.et{1};
    tmptemplate=pm_defs.mflags.template{1};
    pm_defs.mflags.template=tmptemplate;
end
%--------------------------------------------------------------------------
% Load measured field map data - phase and magnitude, real and imaginary or
% precalculated fieldmap
%--------------------------------------------------------------------------
if isfield(job,'precalcfieldmap')
    fm_imgs=spm_vol(job.precalcfieldmap{1});
    if isfield(job,'magfieldmap') && iscell(job.magfieldmap)
        if ~isempty(job.magfieldmap{1})
            pm_defs.magfieldmap=spm_vol(job.magfieldmap{1});
        end
    else
        job.matchvdm=0;
        job.matchanat=0;
        pm_defs.maskbrain=0;
    end
    pm_defs.uflags.iformat='';
elseif isfield(job,'phase') && isfield(job,'magnitude')% && using presub
    tmp=FieldMap('Scale',spm_vol(job.phase{1}));
    fm_imgs=[spm_vol(tmp.fname) spm_vol(job.magnitude{1})];
    pm_defs.uflags.iformat='PM';
elseif isfield(job,'shortphase') && isfield(job,'shortmag')% && using double phase and magnitude
    tmp1=FieldMap('Scale',spm_vol(job.shortphase{1}));
    tmp2=FieldMap('Scale',spm_vol(job.longphase{1}));
    fm_imgs=[spm_vol(tmp1.fname) spm_vol(job.shortmag{1}) spm_vol(tmp2.fname) spm_vol(job.longmag{1})];
    pm_defs.uflags.iformat='PM';
elseif isfield(job,'shortreal') && isfield(job,'shortimag')% && using real & imag
    fm_imgs=[spm_vol(job.shortreal{1}) spm_vol(job.shortimag{1}) spm_vol(job.longreal{1}) spm_vol(job.longimag{1})];
    pm_defs.uflags.iformat='RI';
else
    error('Do not know what to do with this data. Please check your job');
end

%--------------------------------------------------------------------------
% Load epi session data
%--------------------------------------------------------------------------
nsessions=0;
if ~isempty(job.session)
    nsessions=size(job.session,2);
    for sessnum=1:nsessions
        epi_img{sessnum}=job.session(sessnum).epi{1};
    end
else
    epi_img=[];
end

%--------------------------------------------------------------------------
% Load matching, unwarping and session name options
%--------------------------------------------------------------------------
if ~isempty(job.matchvdm)
    pm_defs.match_vdm=job.matchvdm;
else
    pm_defs.match_vdm=0;
end

if ~isempty(job.writeunwarped)
    pm_defs.write_unwarped=job.writeunwarped;
else
    pm_defs.write_unwarped=0;
end

if ~isempty(job.sessname)
    pm_defs.sessname=job.sessname;
else
    pm_defs.sessname='session';
end
%--------------------------------------------------------------------------
% Call FieldMap_create
%--------------------------------------------------------------------------
[VDM, IPcell] = FieldMap_create(fm_imgs,epi_img,pm_defs);

for sessnum=1:max([1 nsessions]);
    
    IP = IPcell{sessnum};
    
    %----------------------------------------------------------------------
    % Display and print results
    %----------------------------------------------------------------------
    fg = spm_figure('FindWin','Graphics');
    if ~isempty(fg)
        spm_figure('Clear','Graphics');
        spm_orthviews('Reset');
    end
    FieldMap('DisplayImage',FieldMap('MakedP'),[.05 .75 .95 .2],1);
    if ~isempty(IP.epiP)
        FieldMap('DisplayImage',IP.epiP,[.05 .5 .95 .2],2);
    end
    if ~isempty(IP.uepiP)
        FieldMap('DisplayImage',IP.uepiP,[.05 .25 .95 .2],3);
    end
    
    %----------------------------------------------------------------------
    % Coregister structural with the unwarped image and display if required
    %----------------------------------------------------------------------
    do_matchanat = 0;
    if iscell(job.anat)
        if ~isempty(job.anat{1})
            IP.nwarp = spm_vol(job.anat{1});
            do_matchanat = job.matchanat;
        end
    end
    
    if ~isempty(IP.nwarp)==1 && ~isempty(IP.epiP)
        if do_matchanat == 1
            fprintf('\nMatching anatomical to unwarped EPI in session %d...\n\n',sessnum);
            FieldMap('MatchStructural',IP);
        end
    end
    
    if ~isempty(IP.nwarp) && ~isempty(IP.epiP)
        FieldMap('DisplayImage',IP.nwarp,[.05 0.0 .95 .2],4);
        % Now need to redisplay other images to make it all look correct
        FieldMap('DisplayImage',FieldMap('MakedP'),[.05 .75 .95 .2],1);
        if ~isempty(IP.epiP)
            FieldMap('DisplayImage',IP.epiP,[.05 .5 .95 .2],2);
        end
        if ~isempty(IP.uepiP)
            FieldMap('DisplayImage',IP.uepiP,[.05 .25 .95 .2],3);
        end
    end
    if ~isempty(fg), spm_print; end
    
    vdm.vdmfile{sessnum} = {VDM{sessnum}.fname};
end
