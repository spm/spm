% This preprocessing script is adapted from 
% script batch_preproc.m by Rik Henson.
% 
% This script will do realign and unwarp with a fieldmap
%
% Chloe Hutton 24/08/05
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Chloe Hutton 
% $Id: Unwarp_batch.m 3692 2010-01-21 21:43:31Z guillaume $

which spm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Enter the name of the directory that contains data/subjects 
%
owd = '/data/FieldMap_examples/sonata_subject';
cd(owd);

defaults = spm('defaults','FMRI');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1) Define subjects and sessions and directory names...
% eg: owd/s1/run1
% 2) Or if images are in current directory, subnames and sessnames 
% can be left empty.
% 3) The vdm_* file should be in a separate directory which is
% within the subject or session directory
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of subjects
nsubs=1;
%subnames={'s1'};
subnames={'images'};

% Sessions for each subject
nsess=[1];
%sessnames={'run1'};
sessnames={''};

P = cell(nsubs,1);

% If including fieldmap pcpm=1;
pcpm=1;

% Name of fieldmap directory
fmapdir='fieldmap'

% These are default flags for coregistration and unwarping which can
% be changed if required

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up flags
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FlagsC  = struct('quality',defaults.realign.estimate.quality,'fwhm',5,'rtm',0,'wrap',[0 0 0],'interp',defaults.realign.estimate.interp);
FlagsR = struct('wrap',[0 0 0],'interp',4,'mask',1,'which',2,'mean',1);

uwe_flags = struct(  'order',       defaults.unwarp.estimate.basfcn,...
            'sfP',         [],...
            'regorder',    defaults.unwarp.estimate.regorder,...
            'lambda',      defaults.unwarp.estimate.regwgt,...
            'jm',          defaults.unwarp.estimate.jm,...
            'fot',         [4 5],...    % pitch and roll
            'sot',         [],...       % no second-order 
            'fwhm',        defaults.unwarp.estimate.fwhm,...
            'rem',         defaults.unwarp.estimate.rem,...
            'noi',         defaults.unwarp.estimate.noi,...
            'exp_round',   defaults.unwarp.estimate.expround);
        
uwr_flags = struct(  'wrap',        [0 0 0],...
            'interp',      4,...
            'mask',        1,...
            'which',       2,...
            'mean',        1,...
            'udc',         1);
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:nsubs
      
   pp = cell(1,nsess(i));

   if pcpm;ppm = cell(1,nsess(i)); else ppm = []; end % Are fieldmaps included

   for ses=1:nsess(i)
        dir{ses} = sprintf('%s%s%s%s%s',owd,filesep,subnames{i},filesep,sessnames{ses});
        pp{ses}= spm_get('Files',dir{ses},'f*.img');

        % Get fieldmaps from fmapdir
        % NB the naming here will depend whether 
        % fieldmap is session specific or not.

        fdir{ses}= deblank(sprintf('%s%s%s%s%s%s%s',owd,filesep,subnames{i},filesep,sessnames{ses},filesep,fmapdir));
        pm = spm_get('Files',fdir{ses},'vdm_*.img');
        if pcpm & ~isempty(pm)
           ppm{ses} = pm;
    elseif pcpm & isempty(pm) & ~isempty(ppm{ses-1})
       ppm{ses} = ppm{ses-1};
    elseif pcpm & isempty(pm) & ses==1
           msg=sprintf('Cannot find a vdm file for %s%s%s',subnames{i},filesep,sessnames{ses});
           error(msg);
        end
   end
   P{i}=pp;
   pP{i} = ppm;
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Realignment 
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   disp(sprintf('sub %d, coregistering...',i))
    
   spm_realign(P{i},FlagsC);
            
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Unwarping
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
   disp(sprintf('sub %d, unwarping...',i))
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % Get the space of the first image of first session for each subject
   clear ads
   tmpP = spm_vol(P{i}{1}(1,:));
   uwe_flags.M  = tmpP.mat;    

   for j=1:length(P{i})
      
      % Are we including a field map?
      if pcpm > 0
         uwe_flags.sfP = pP{i}{j};
      end       
      ds = spm_uw_estimate(P{i}{j},uwe_flags);
      ads(j) = ds;
      [path,name,ext,ver] = fileparts(P{i}{j}(1,:));
      pefile = fullfile(path,[name '_uw.mat']);
      save(pefile,'ds');
   end
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
   spm_uw_apply(ads,uwr_flags);
        
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
end
