function VDM=FieldMap_preprocess(fm_dir,epi_dir,pm_defs)
%
% Function to prepare fieldmap data for processing
% 
% FORMAT VDM = FieldMap_preprocess(fm_dir,epi_dir,pm_defs)
% fm_dir    - name of directory containing fieldmap images
% epi_dir   - name of directory containing epi images (needs first epi in time
%             series to match the fieldmap to).
%
% pm_defs   - vector containing following values (optional flags in brackets): 
%             [te1,te2,epifm,tert,kdir,(mask),(match)];
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
% VDM       - file pointer to the VDM file (voxel displacement map) required
%             for the Unwarping process. This will be written to the
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
% VDM=FieldMap_preprocess(fm_dir,epi_dir,[10.0,14.76,0,32,1]);
%
% Sonata EPI fieldmap parameters and default EPI fMRI
% VDM=FieldMap_preprocess(fm_dir,epi_dir,[25,34.5,1,32,1]);
%
% Allegra Siemens fieldmap parameters and default EPI fMRI
% VDM=FieldMap_preprocess(fm_dir,epi_dir,[10.0,12.46,0,21,1]);
%
% Allegra EPI fieldmap parameters and default EPI fMRI
% VDM=FieldMap_preprocess(fm_dir,epi_dir,[19,29,1,21,1]);
%  
% It is also possible to switch off the brain masking which is
% done by default with a siemens field map (set 6th flag to 0) 
% and the matching of the fieldmap to the EPI (set 7th flag to 0).
% 
%_________________________________________________________________
% FieldMap_preprocess.m                           Chloe Hutton 23/08/05
%
if nargin < 3
  error('Usage: FieldMap_preprocess(fm_dir,epi_dir,pm_defs)');
end

if size(pm_defs,1)<5 & size(pm_defs,2)<5
  error('The following parameters are required: te1, te2, epifm, tert, kdir');
end

pm_defaults;

% Update default parameters
pm_def.SHORT_ECHO_TIME= pm_defs(1);
pm_def.LONG_ECHO_TIME=pm_defs(2);
pm_def.EPI_BASED_FIELDMAPS=pm_defs(3);
pm_def.TOTAL_EPI_READOUT_TIME=pm_defs(4);
pm_def.K_SPACE_TRAVERSAL_BLIP_DIR=pm_defs(5);

% If using a non-epi fieldmap, the input format will be 'PM'
if pm_def.EPI_BASED_FIELDMAPS==1
   pm_def.INPUT_DATA_FORMAT='RI';
   pm_def.MASK=0;

% Do brain masking for Siemens fieldmaps unless switched off in pm_defs(6)
elseif pm_def.EPI_BASED_FIELDMAPS==0
   pm_def.INPUT_DATA_FORMAT='PM';
   pm_def.MASK=1;
   if size(pm_defs,1)>5 | size(pm_defs,2)>5
      if pm_defs(6)==1
         pm_def.MASK=1;
      elseif pm_defs(6)==0
         pm_def.MASK=0;
      end
   end
else
   error('Sorry the parameter epifm must be 0 or 1');
end

% Match the VDM to EPI unless unless switched off in pm_defs(7)
pm_def.match_vdm=1;
if size(pm_defs,1)>6 | size(pm_defs,2)>6
   if pm_defs(7)==0
      pm_def.match_vdm=0;
   end
end


%----------------------------------------------------------------------
% Load epi data from data directory
%----------------------------------------------------------------------
epi_img = spm_get('Files',epi_dir,'f*.img');
epi_img=epi_img(1,:);

%----------------------------------------------------------------------
% Load field map data from fieldmap directory
%----------------------------------------------------------------------
if pm_def.INPUT_DATA_FORMAT=='PM'
   fm_imgs = spm_get('Files',fm_dir,'s*.img');
   if ~isempty(fm_imgs)
      nfiles=size(fm_imgs,1);
      if nfiles~=3
         msg=sprintf('Wrong number of field map (s*.img) images! There should be 3!');
         error(msg);
      else
         % Need first of two mag images and the phase
         % These may have been acquired in either order so need to check for this
         nn=findstr(fm_imgs(1,:),'-');
         if isempty(findstr(fm_imgs(2,:),fm_imgs(1,nn(1):nn(2))))
            phase=fm_imgs(1,:);
            mag=fm_imgs(2,:);
         else
            mag=fm_imgs(1,:);
            phase=fm_imgs(3,:);
         end
         scphase=FieldMap('Scale',phase);
         fm_imgs=str2mat(scphase.fname,mag);
      end
   else
      msg=sprintf('Sorry. I cannot find the fieldmap files in %s',fm_dir);
      error(msg);
   end       
elseif pm_def.INPUT_DATA_FORMAT=='RI'  

   % This expects to find six EPI field map files: 
   % 3 short (real, imag and mag) and 3 long (real, imag and mag).

   all_fm_imgs = spm_get('Files',fm_dir,'f*.img');
   nfiles=size(all_fm_imgs,1);
   if nfiles~=6
         msg=sprintf('Wrong number of field map (f*.img) images! There should be 6!');
         error(msg);
   else

      % Now the FieldMap creation expects the files in the order:
      % short-real, short-imag, long-real, long-imag

      fm_imgs=all_fm_imgs([2 1 5 4],:); 
      if (isempty(strfind(fm_imgs(1,:),'short-real')) | isempty(strfind(fm_imgs(2,:),'short-imag')) | isempty(strfind(fm_imgs(3,:),'long-real')) | isempty(strfind(fm_imgs(4,:),'long-imag')))
         msg=sprintf('There is a problem with the files. FieldMap needs short-real, short-imag, long-real, long-imag');
         error(msg);
      end
   end
else
   error('Funny input format specified. FieldMap needs short-real, short-imag, long-real, long-imag')
end

%----------------------------------------------------------------------
% Run function to create vdm file 
%----------------------------------------------------------------------
VDM = FieldMap_create(fm_imgs,epi_img,pm_def);

return
