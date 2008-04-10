function vdm=Fieldmap_Run(job)
% Auxillary file for running FieldMap jobs
%
% FORMAT vdm = Fieldmap_Run(job)
%
% job  - FieldMap job structure containing following elements
%        common to all jobs and specific to the type of job
%     Common to all jobs:
%        defaults - cell array containing name string of the defaults file
%        options - structure containing the following:
%           display - display results or not (1/0)
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
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Chloe Hutton & Jesper Andersson
% $Id: FieldMap_Run.m 1358 2008-04-10 11:20:26Z guillaume $
%_________________________________________________________________

%
%----------------------------------------------------------------------  
% Set up default parameters and structures 
%----------------------------------------------------------------------
spm('defaults','FMRI');

% Open the FieldMap control window with visibility off. This allows the
% graphics display to work.
FieldMap('Welcome','Off');
IP = FieldMap('Initialise'); % Gets default params from pm_defaults

% Here load the selected defaults file necessary
m_file = job.defaults;
m_file = m_file{1}(1:end-2);
m_file = spm_str_manip(m_file,'t');
IP = FieldMap('SetParams',m_file); % Gets default params from pm_defaults

% Load precalculated Hz fieldmap if required
if isfield(job,'precalcfieldmap')
    IP.pP = spm_vol(job.precalcfieldmap{1});
end
%----------------------------------------------------------------------
% Load measured field map data - phase and magnitude or real and imaginary
%----------------------------------------------------------------------
if ~isempty(IP.pP)
   IP.fm.fpm = spm_read_vols(IP.pP);
   IP.fm.jac = pm_diff(IP.fm.fpm,2);
   if isfield(IP,'P') & ~isempty(IP.P{1})
      IP.P = cell(1,4);
   end
elseif IP.uflags.iformat=='PM' &  isfield(job,'phase') % && using presub 
   IP.P{1}=spm_vol(job.phase{1});
   tmp=FieldMap('Scale',IP.P{1}.fname);
   IP.P{1} = spm_vol(tmp.fname);
   IP.P{2}=spm_vol(job.magnitude{1});
elseif IP.uflags.iformat=='PM' &  isfield(job,'shortphase') % && using double phase and magnitude
   IP.P{1}=spm_vol(job.shortphase{1});
   tmp=FieldMap('Scale',IP.P{1}.fname);
   IP.P{2}=spm_vol(job.shortmag{1});
   IP.P{3}=spm_vol(job.longphase{1});
   tmp=FieldMap('Scale',IP.P{3}.fname);
   IP.P{4}=spm_vol(job.longmag{1});
elseif IP.uflags.iformat=='RI' 
   IP.P{1}=spm_vol(job.shortreal{1});
   IP.P{2}=spm_vol(job.shortimag{1});
   IP.P{3}=spm_vol(job.longreal{1});
   IP.P{4}=spm_vol(job.longimag{1});
else
    error('Do not know what to do with this data. Please check your defaults file');
end

%----------------------------------------------------------------------
% Load job options
%----------------------------------------------------------------------

% Only do matchvdm and writeunwarped if an epi has been selected
do_matchvdm = 0;
do_write_unwarped = 0;
if ~isempty(job.options.epi)      
    IP.epiP = spm_vol(job.options.epi{1}); 
    do_matchvdm = job.options.matchvdm;
    do_writeunwarped = job.options.writeunwarped;
end

% Only do matchanat if anatomical image has been selected
do_matchanat = 0;
if ~isempty(job.options.anat{1})      
    IP.nwarp = spm_vol(job.options.anat{1});
    do_matchanat = job.options.matchanat;
end

do_display = job.options.display;

%----------------------------------------------------------------------
% Create field map (in Hz) - this routine calls the unwrapping
%----------------------------------------------------------------------
if isempty(IP.pP)

   IP.fm = FieldMap('CreateFieldMap',IP);

%----------------------------------------------------------------------
% Write out field map
% Outputs -> fpm_NAME-OF-FIRST-INPUT-IMAGE.img
%----------------------------------------------------------------------

   FieldMap('Write',IP.P{1},IP.fm.fpm,'fpm_',64,'Fitted phase map in Hz');
end

%----------------------------------------------------------------------
% Convert Hz to voxels and write voxel displacement map 
% Outputs -> vdm5_NAME-OF-FIRST-INPUT-IMAGE.img
%----------------------------------------------------------------------

[IP.vdm, IP.vdmP]=FieldMap('FM2VDM',IP);

%----------------------------------------------------------------------
% Match VDM and/or unwarp EPI 
%----------------------------------------------------------------------
if ~isempty(IP.epiP) 
      
    if do_matchvdm == 1

%----------------------------------------------------------------------
% Match voxel displacement map to image
% Outputs -> mag_NAME-OF-FIRST-INPUT-IMAGE.img
%----------------------------------------------------------------------
        IP.vdmP = FieldMap('MatchVDM',IP);
    end    
    
%----------------------------------------------------------------------
% Unwarp the EPI
% Outputs -> mag_NAME-OF-FIRST-INPUT-IMAGE.img
%----------------------------------------------------------------------    

    IP.uepiP = FieldMap('UnwarpEPI',IP);
    
    if do_writeunwarped == 1
        
%----------------------------------------------------------------------
% Outputs -> uNAME-OF-EPI.img
%----------------------------------------------------------------------
       unwarp_info=sprintf('Unwarped EPI:echo time difference=%2.2fms, EPI readout time=%2.2fms, Jacobian=%d',IP.uflags.etd, IP.tert,IP.ajm);    
       IP.uepiP = FieldMap('Write',IP.epiP,IP.uepiP.dat,'u',IP.epiP.dt(1),unwarp_info);

    end
end

if ~isempty(IP.nwarp) == 1

    if do_matchanat == 1
        
%----------------------------------------------------------------------
% Coregister structural with the unwarped image
%----------------------------------------------------------------------
       FieldMap('MatchStructural',IP);
    end
end

if do_display == 1

    FieldMap('DisplayImage',FieldMap('MakedP'),[.05 .75 .95 .2],1);
    if ~isempty(IP.epiP)
       FieldMap('DisplayImage',IP.epiP,[.05 .5 .95 .2],2);
    end
    if ~isempty(IP.uepiP)
       FieldMap('DisplayImage',IP.uepiP,[.05 .25 .95 .2],3);
    end
    if ~isempty(IP.nwarp)
       FieldMap('DisplayImage',IP.nwarp,[.05 0.0 .95 .2],4);
    end
end
        
%______________________________________________________________________
vdm=IP.vdmP;
