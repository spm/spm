function VDM = FieldMap_preprocess(fm_dir,epi_dir,pm_defs,sessname)
% Function to prepare fieldmap data for processing
% 
% FORMAT VDM = FieldMap_preprocess(fm_dir,epi_dir,pm_defs,sessname)
% fm_dir    - name of directory containing fieldmap images
% epi_dir   - name of directory containing epi images (needs first epi in time
%             series to match the fieldmap to).
%             This can also be a cell array of directory names for creating
%             session-specific versions of a vdm file where each vdm file
%             is matched to the first image of each EPI directory specified.
%             Each session specific vdm file will be saved with the name
%             vdm5_XXXX_'sessname'N.img where 'sessname is 'session' by
%             default or another name if specified by the user as the fourth
%             argument to the script.
% pm_defs   - vector containing following values (optional flags in brackets): 
%             [te1,te2,epifm,tert,kdir,(mask),(match),(write)];
%
% te1       - short echo time
% te2       - long echo time
% epifm     - epi-based fieldmap (1/0)?
% tert      - total echo readout time
% kdir      - blip direction (+1/-1)
% mask      - (optional flag, default=1) Do brain masking or not 
%             (only if non-epi fieldmap)
% match     - (optional flag, default=1) Match fieldmap to epi or not
%
% writeunwarped 
%           - (optional flag, default=1) Write unwarped epi or not
%
% sessname  - (optional string, default='session') This will be the name
%             extension followed by an incremented integer for session specific vdm files.
%
% VDM       - cell array of file pointers to the VDM file(s) (voxel displacement map) 
%             required for the Unwarping process. This will be written to the
%             same directory as the fieldmap data.
%
% NB:
% 1) This function takes input directory names and parameters and puts them 
% into the correct format for creating fieldmaps
% 2) The function assumes that only the fieldmap images are in the
% fieldmap directory
% 
% Below is a list of the most common sequences and parameter lists 
% used at the FIL:
%
% Sonata Siemens fieldmap parameters and default EPI fMRI'); 
% VDM=FieldMap_preprocess(fm_dir,epi_dir,[10.0,14.76,0,32,-1]);
%
% Allegra Siemens fieldmap parameters and default EPI fMRI
% VDM=FieldMap_preprocess(fm_dir,epi_dir,[10.0,12.46,0,21.12,-1]);
%
% Allegra Siemens fieldmap parameters and extended FOV EPI fMRI
% VDM=FieldMap_preprocess(fm_dir,epi_dir,[10.0,12.46,0,23.76,-1]);
%
% Allegra Siemens fieldmap parameters and 128 EPI fMRI
% VDM=FieldMap_preprocess(fm_dir,epi_dir,[10.0,12.46,0,71.68,-1]);
%  
% It is also possible to switch off the brain masking which is
% done by default with a siemens field map (set 6th flag to 0) 
% and the matching of the fieldmap to the EPI (set 7th flag to 0).
% 
% This function generates session specific versions of the vdm file that
% have been matched to the first image of each session.
%__________________________________________________________________________
% Copyright (C) 2006-2014 Wellcome Trust Centre for Neuroimaging

% Chloe Hutton 
% $Id: FieldMap_preprocess.m 5962 2014-04-17 12:47:43Z spm $

 
if nargin < 3
  error('Usage: FieldMap_preprocess(fm_dir,epi_dir,pm_defs)');
end

if size(pm_defs,1)<5 && size(pm_defs,2)<5
  error('The following parameters are required: te1, te2, epifm, tert, kdir');
end

pm_defaults;
if nargin<4
  pm_def.sessname='session';
else
  pm_def.sessname=sessname;
end

% Update default parameters
pm_def.SHORT_ECHO_TIME= pm_defs(1);
pm_def.LONG_ECHO_TIME=pm_defs(2);
pm_def.EPI_BASED_FIELDMAPS=pm_defs(3);
pm_def.TOTAL_EPI_READOUT_TIME=pm_defs(4);
pm_def.K_SPACE_TRAVERSAL_BLIP_DIR=pm_defs(5);

% If using a non-epi fieldmap, the input format will be 'PM'
if pm_def.EPI_BASED_FIELDMAPS==1
   pm_def.INPUT_DATA_FORMAT='RI';
   pm_def.MASKBRAIN=0;

% Do brain masking for Siemens fieldmaps unless switched off in pm_defs(6)
elseif pm_def.EPI_BASED_FIELDMAPS==0
   pm_def.INPUT_DATA_FORMAT='PM';
   pm_def.MASKBRAIN=1;
   if size(pm_defs,1)>5 || size(pm_defs,2)>5
      if pm_defs(6)==1
         pm_def.MASKBRAIN=1;
      elseif pm_defs(6)==0
         pm_def.MASKBRAIN=0;
      end
   end
else
   error('Sorry the parameter epifm must be 0 or 1');
end

%--------------------------------------------------------------------------
% Match the VDM to EPI unless unless switched off in pm_defs(7)
%--------------------------------------------------------------------------
pm_def.match_vdm=1;
if size(pm_defs,1)>6 || size(pm_defs,2)>6
   if pm_defs(7)==0
      pm_def.match_vdm=0;
   end
end

%--------------------------------------------------------------------------
% Write the unwarped EPI unless unless switched off in pm_defs(8)
%--------------------------------------------------------------------------
pm_def.write_unwarped=1;
if size(pm_defs,1)>7 || size(pm_defs,2)>7
   if pm_defs(8)==0
      pm_def.write_unwarped=0;
   end
end

%--------------------------------------------------------------------------
% Load epi data from data directory for each session if necessary
%--------------------------------------------------------------------------
if iscell(epi_dir)
   nsessions = size(epi_dir,2);
   for sessnum=1:nsessions
      epi_all = strvcat(spm_select('List',epi_dir{sessnum},'^f.*\.nii$'),spm_select('List',epi_dir{sessnum},'^f.*\.img$'));
      epi_img{sessnum}=fullfile(epi_dir{sessnum},epi_all(1,:));
   end
elseif ischar(epi_dir)
   nsessions=1;
   sessnum=1;
   epi_all = strvcat(spm_select('List',epi_dir,'^f.*\.nii$'),spm_select('List',epi_dir,'^f.*\.img$'));
   epi_img{sessnum}=fullfile(epi_dir,epi_all(1,:));
end

%--------------------------------------------------------------------------
% Load field map data from fieldmap directory
%--------------------------------------------------------------------------
if strcmp(pm_def.INPUT_DATA_FORMAT,'PM')
   fm_imgs = strvcat(spm_select('List',fm_dir,'^s.*\.nii$'),spm_select('List',fm_dir,'^s.*\.img$'));
   if ~isempty(fm_imgs)
      nfiles=size(fm_imgs,1);
      % Added the next few lines so the scaled fieldmap files aren't picked
      % up
      if nfiles>3
          nnfiles=0;
          for filenum=1:nfiles        
              if ~strncmp(deblank(fm_imgs(filenum,:)),'sc',2)
                  nfm_imgs(nnfiles+1,:)=fm_imgs(filenum,:);
                  nnfiles=nnfiles+1;
              end
          end
          nfiles=nnfiles;
          fm_imgs=nfm_imgs;
      end
      if nfiles~=3
         error('Wrong number of field map (s*.img) images! There should be 3!');
      else
         % Need first of two mag images and the phase
         % These may have been acquired in either order so need to check for this
         nn=findstr(fm_imgs(1,:),'-');
         if isempty(findstr(fm_imgs(2,:),fm_imgs(1,nn(1):nn(2))))
            phase=fullfile(fm_dir,fm_imgs(1,:));
            mag=fullfile(fm_dir,fm_imgs(2,:));
         else
            mag=fullfile(fm_dir,fm_imgs(1,:));
            phase=fullfile(fm_dir,fm_imgs(3,:));
         end
         scphase=FieldMap('Scale',phase);
         fm_imgs=char(scphase.fname,mag);
      end
   else
      error('Sorry. I cannot find the fieldmap files in %s',fm_dir);
   end       
elseif strcmp(pm_def.INPUT_DATA_FORMAT,'RI')

   % This expects to find six EPI field map files: 
   % 3 short (real, imag and mag) and 3 long (real, imag and mag).

   all_fm_imgs = strvcat(spm_select('List',fm_dir,'^f.*\.nii$'), spm_select('List',fm_dir,'^f.*\.img$'));
   nfiles=size(all_fm_imgs,1);
   if nfiles~=6
         error('Wrong number of field map (f*.img) images! There should be 6!');
   else

      % Now the FieldMap creation expects the files in the order:
      % short-real, short-imag, long-real, long-imag

      fm_imgs=all_fm_imgs([2 1 5 4],:); 
      if (isempty(strfind(fm_imgs(1,:),'short-real')) || isempty(strfind(fm_imgs(2,:),'short-imag')) || isempty(strfind(fm_imgs(3,:),'long-real')) | isempty(strfind(fm_imgs(4,:),'long-imag')))
         error('There is a problem with the files. FieldMap needs short-real, short-imag, long-real, long-imag');
      end
   end
else
   error('Funny input format specified. FieldMap needs short-real, short-imag, long-real, long-imag')
end

%--------------------------------------------------------------------------
% Run function to create vdm file 
%--------------------------------------------------------------------------
VDM = FieldMap_create(fm_imgs,epi_img,pm_def);
