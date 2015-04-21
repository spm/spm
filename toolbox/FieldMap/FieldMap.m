function varargout = FieldMap(varargin)
% FieldMap is an SPM Toolbox for creating field maps and unwarping EPI.
% A full description of the toolbox and a usage manual can be found in
% FieldMap.man. This can launched by the toolbox help button or using
% `spm_help FieldMap.man`. The theoretical and practical principles behind
% the toolbox are described in FieldMap_principles.man.
%
% FORMAT FieldMap
%
% FieldMap launches the GUI-based toolbox. Help is available via the help
% button (which calls spm_help FieldMap.man). FieldMap is a multi function
% function so that the toolbox routines can also be accessed without using
% the GUI. A description of how to do this can be found in FieldMap_ngui.m
% 
% Input parameters and the mode in which the toolbox works can be 
% customised using the defaults file called pm_defaults.m. 
% 
% Main data structure:
%
% IP.P              : 4x1 cell array containing real part short TE,
%                     imaginary part short TE, real part long TE and
%                     imaginary part long TE.
% IP.pP             : Cell containing pre-calculated phase map. N.B.
%                     IP.P and IP.pP are mutually exclusive.
% IP.epiP           : Cell containing EPI image used to demonstrate 
%                     effects of unwarping.
% IP.fmagP          : Cell containing fieldmap magnitude image used for 
%                     coregistration
% IP.wfmagP         : Cell containing forward warped fieldmap magnitude 
%                     image used for coregistration
% IP.uepiP          : Cell containing unwarped EPI image.
% IP.nwarp          : Cell containing non-distorted image.
% IP.vdmP           : Cell containing the voxel displacement map (VDM)
% IP.et             : 2x1 Cell array with short and long echo-time (ms).
% IP.epifm          : Flag indicating EPI based field map (1) or not (0).
% IP.blipdir        : Direction of phase-encode blips for k-space traversal
%                    (1 = positive or -1 = negative)
% IP.ajm            : Flag indicating if Jacobian modulation should be applied 
%                    (1) or not (0).
% IP.tert           : Total epi readout time (ms).
% IP.maskbrain      : Flag indicating whether to mask the brain for fieldmap creation
% IP.uflags         : Struct containing parameters guiding the unwrapping.
%                     Further explanations of these parameters are in 
%                     FieldMap.man and pm_make_fieldmap.m 
% .iformat          : 'RI' or 'PM' 
% .method           : 'Huttonish', 'Mark3D' or 'Mark2D'
% .fwhm             : FWHM (mm) of Gaussian filter for field map smoothing
% .pad              : Size (in-plane voxels) of padding kernel. 
% .etd              : Echo time difference (ms).   
% .bmask              
%
% IP.mflags         : Struct containing parameters for brain maskin
% .template         : Name of template for segmentation.
% .fwhm             : fwhm of smoothing kernel for generating mask.
% .nerode           : number of erosions
% .ndilate          : number of dilations
% .thresh           : threshold for smoothed mask. 
% .reg              : bias field regularisation
% .graphics         : display or not
%
% IP.fm             : Struct containing field map information
% IP.fm.upm         : Phase-unwrapped field map (Hz).
% IP.fm.mask        : Binary mask excluding the noise in the phase map.
% IP.fm.opm         : "Raw" field map (Hz) (not unwrapped).
% IP.fm.fpm         : Phase-unwrapped, regularised field map (Hz).
% IP.fm.jac         : Partial derivative of field map in y-direction.
%
% IP.vdm            : Struct with voxel displacement map information
% IP.vdm.vdm        : Voxel displacement map (scaled version of IP.fm.fpm).
% IP.vdm.jac        : Jacobian-1 of forward transform.
% IP.vdm.ivdm       : Inverse transform of voxel displacement 
%                     (used to unwarp EPI image if field map is EPI based)
%                     (used to warp flash image prior to coregistration
%                     when field map is flash based (or other T2 weighting).
% IP.vdm.ijac       : Jacobian-1 of inverse transform.
% IP.jim            : Jacobian sampled in space of EPI.
%
% IP.cflags         : Struct containing flags for coregistration 
%                     (these are the default SPM coregistration flags - 
%                     defaults.coreg).
% .cost_fun         
% .sep             
% .tol
% .fwhm    
%
%__________________________________________________________________________
% Refs and Background reading:
%
% Jezzard P & Balaban RS. 1995. Correction for geometric distortion in
% echo planar images from Bo field variations. MRM 34:65-73.
%
% Hutton C et al. 2002. Image Distortion Correction in fMRI: A Quantitative 
% Evaluation, NeuroImage 16:217-240.
%
% Cusack R & Papadakis N. 2002. New robust 3-D phase unwrapping
% algorithms: Application to magnetic field mapping and
% undistorting echoplanar images. NeuroImage 16:754-764.
%
% Jenkinson M. 2003. Fast, automated, N-dimensional phase-
% unwrapping algorithm. MRM 49:193-197.
%__________________________________________________________________________
% Acknowledegments:
% 
% Wellcome Trust and IBIM Consortium
%__________________________________________________________________________
% Copyright (C) 2006-2014 Wellcome Trust Centre for Neuroimaging

% Jesper Andersson and Chloe Hutton 
% $Id: FieldMap.m 6416 2015-04-21 15:34:10Z guillaume $


persistent PF FS WS PM   % GUI related constants
persistent ID            % Image display
persistent IP            % Input and results
persistent DGW           % Delete Graphics Window (if we created it)
global st                % Global for spm_orthviews

col = [0.9 0.9 0.9];

if nargin == 0
   Action = 'welcome';
else
   Action = varargin{1};
end

% Actions are divided into 3 categories:
% 1) Functions that can be called directly from the GUI are at the
% beginning. 
% 2) Next come other GUI-related 'utility' functions - ie those that
% care of the GUI behind the scenes.
% 3) Finally, call-back functions that are not GUI dependent and can 
% be called from a script are situated towards the end. 
% See FieldMap_ngui.m for details on how to use these.

switch lower(Action)
%=======================================================================
%
% Create and initialise GUI for FieldMap toolbox
%
%=======================================================================

   case 'welcome'
              
      % Unless specified, set visibility to on
      if nargin==2
          if ~strcmp(varargin{2},'off') && ~strcmp(varargin{2},'Off')
            visibility = 'On';
          else
            visibility = 'Off';
          end
      else
          visibility = 'On';
      end
      
      DGW = 0;

      % Get all default values (these may effect GUI)
      IP = FieldMap('Initialise');
      %
      % Since we are using persistent variables we better make sure
      % there is no-one else out there.
      %
      if ~isempty(PM)
         figure(PM);
         set(PM,'Visible',visibility);
         return
      end
 
      S0   = spm('WinSize','0',1);
      WS   = spm('WinScale');   
      FS   = spm('FontSizes');
      PF   = spm_platform('fonts');
      rect = [100 100 510 520].*WS;

      PM = figure('IntegerHandle','off',...
       'Name','FieldMap Toolbox (version 2.1)',...
       'NumberTitle','off',...
       'Tag','FieldMap',...
       'Position',[S0(1),S0(2),0,0] + rect,...
       'Resize','off',...
       'Pointer','Arrow',...
       'Color',[1 1 1]*.8,...
       'MenuBar','none',...
       'DefaultTextFontName',PF.helvetica,...
       'DefaultTextFontSize',FS(10),...
       'DefaultUicontrolFontName',PF.helvetica,...
       'DefaultUicontrolFontSize',FS(10),...
       'HandleVisibility','on',...
       'Visible',visibility,...
       'DeleteFcn','FieldMap(''Quit'');');

%-Create phase map controls

%-Frames and text

      uicontrol(PM,'Style','Frame','BackgroundColor',col,...
                'Position',[10 270 490 240].*WS);
      uicontrol(PM,'Style','Text','String','Create field map in Hz',...
            'Position',[100 480 300 020].*WS,...
            'ForegroundColor','k','BackgroundColor',col,...
                'FontName',PF.times,'FontAngle','Italic')
      uicontrol(PM,'Style','Text','String','Short TE',...
            'Position',[25 452 60 20].*WS,...
            'ForegroundColor','k','BackgroundColor',col,...
                'FontName',PF.times)
      uicontrol(PM,'Style','Text','String','ms',...
            'Position',[78 430 30 20].*WS,...
            'ForegroundColor','k','BackgroundColor',col,...
                'FontName',PF.times)
      uicontrol(PM,'Style','Text','String','Long TE',...
            'Position',[25 392 60 20].*WS,...
            'ForegroundColor','k','BackgroundColor',col,...
                'FontName',PF.times)
      uicontrol(PM,'Style','Text','String','ms',...
            'Position',[78 370 30 20].*WS,...
            'ForegroundColor','k','BackgroundColor',col,...
                'FontName',PF.times)
      uicontrol(PM,'Style','Text','String',{'Mask brain:'},...
            'Position',[30 310 80 35].*WS,...
            'ForegroundColor','k','BackgroundColor',col,...
                'FontName',PF.times)
      uicontrol(PM,'Style','Text','String',{'Precalculated','field map:'},...
            'Position',[30 275 80 35].*WS,...
            'ForegroundColor','k','BackgroundColor',col,...
                'FontName',PF.times)
      uicontrol(PM,'Style','Text','String',{'Field map value:'},...
            'Position',[240 280 100 20].*WS,...
            'ForegroundColor','k','BackgroundColor',col,...
                'FontName',PF.times)
      uicontrol(PM,'Style','Text','String',{'Hz'},...
            'Position',[403 280 20 20].*WS,...
            'ForegroundColor','k','BackgroundColor',col,...
                'FontName',PF.times)


%-Objects with Callbacks

      uicontrol(PM,'Style','radiobutton',...
          'String','RI',...
          'Position',[30 480 44 022].*WS,...
          'ToolTipString','Select Real/Imaginary input images',...
          'CallBack','FieldMap(''InputFormat'');',...
          'UserData',1,...
          'BackgroundColor',col,...
          'Tag','IFormat');
      uicontrol(PM,'Style','radiobutton',...
          'String','PM',...
          'Position',[76 480 44 022].*WS,...
          'ToolTipString','Select phase/magnitude input images',...
          'CallBack','FieldMap(''InputFormat'');',...
          'UserData',2,...
          'BackgroundColor',col,...
          'Tag','IFormat');
      uicontrol(PM,'Style','Edit',...
          'Position',[30 430 50 024].*WS,...
          'ToolTipString','Give shortest echo-time in ms',...
          'CallBack','FieldMap(''EchoTime'');',...
          'UserData',1,...
          'Tag','ShortTime');
      uicontrol(PM,'Style','Edit',...
          'Position',[30 370 50 024].*WS,...
          'ToolTipString','Give longest echo-time in ms',...
          'CallBack','FieldMap(''EchoTime'');',...
          'UserData',2,...
          'Tag','LongTime');

      % Here we set up appropriate gui

      FieldMap('IFormat_Gui');

      uicontrol(PM,'String','Calculate',...
          'Position',[265 337 90 022].*WS,...
          'Enable','Off',...
          'ToolTipString','Go ahead and create unwrapped field map (in Hz)',...
          'CallBack','FieldMap(''CreateFieldmap_Gui'');',...
          'Tag','CreateFieldMap');
      uicontrol(PM,'String','Write',...
          'Position',[370 337 90 022].*WS,...
          'Enable','Off',...
          'ToolTipString','Write Analyze image containing fininshed field map (in Hz)',...
          'CallBack','FieldMap(''WriteFieldMap_Gui'');',...
          'Tag','WriteFieldMap');
      uicontrol(PM,'String','Load',...
          'Position',[125 280 90 022].*WS,...
          'ToolTipString','Load Analyze image containing finished field map (in Hz)',...
          'CallBack','FieldMap(''LoadFieldMap_Gui'');',...
          'Tag','LoadFieldMap');  
   
      % Empty string written to Fieldmap value box
      uicontrol(PM,'Style','Text',...
          'Position',[340 280 50 020].*WS,...
          'HorizontalAlignment','left',...
          'String','');

%-Create voxel-displacement-map controls

%-Frames and text

      uicontrol(PM,'Style','Frame','BackgroundColor',col,...
          'Position',[10 10 490 250].*WS);
      uicontrol(PM,'Style','Text',...
          'String','Create voxel displacement map (VDM) and unwarp EPI',...
          'Position',[70 235 350 020].*WS,...
          'ForegroundColor','k','BackgroundColor',col,...
          'FontName',PF.times,'FontAngle','Italic')
      uicontrol(PM,'Style','Text',...
          'String','EPI-based field map',...
          'Position',[50 200 200 020].*WS,...
          'ForegroundColor','k','BackgroundColor',col,...
          'FontName',PF.times)
      uicontrol(PM,'Style','Text',...
          'String','Polarity of phase-encode blips',...
          'Position',[50 170 200 020].*WS,...
          'ForegroundColor','k','BackgroundColor',col,...
          'FontName',PF.times)
      uicontrol(PM,'Style','Text',...
          'String','Apply Jacobian modulation',...
          'Position',[50 140 200 020].*WS,...
          'ForegroundColor','k','BackgroundColor',col,...
          'FontName',PF.times)
      uicontrol(PM,'Style','Text',...
          'String','Total EPI readout time',...
          'Position',[50 110 200 020].*WS,...
          'ForegroundColor','k','BackgroundColor',col,...
          'FontName',PF.times)
      uicontrol(PM,'Style','Text',...
          'String','ms',...
          'Position',[370 110 30 020].*WS,...
          'ForegroundColor','k','BackgroundColor',col,...
          'FontName',PF.times)

%-Objects with Callbacks
       

      uicontrol(PM,'style','RadioButton',...
          'String','Yes',...
          'Position',[300 200 60 022].*WS,...
          'ToolTipString','Field map is based on EPI data and needs inverting before use',...
          'CallBack','FieldMap(''RadioButton'');',...
          'Tag','EPI',...
          'UserData',1);
      uicontrol(PM,'style','RadioButton',...
          'String','No',...
          'Position',[370 200 60 022].*WS,...
          'ToolTipString','Phase-map is based on non-EPI (Flash) data and can be used directly',...
          'CallBack','FieldMap(''RadioButton'');',...
          'Tag','EPI',...
          'UserData',2);

      uicontrol(PM,'style','RadioButton',...
          'String','+ve',...
          'Position',[370 170 60 022].*WS,...
          'ToolTipString','K-space is traversed using positive phase-encode blips',...
          'CallBack','FieldMap(''RadioButton'');',...
          'Tag','BlipDir',...
          'UserData',1);
      uicontrol(PM,'style','RadioButton',...
          'String','-ve',...
          'Position',[300 170 60 022].*WS,...
          'ToolTipString','K-space is traversed using negative phase-encode blips',...
          'CallBack','FieldMap(''RadioButton'');',...
          'Tag','BlipDir',...
          'UserData',2);
      
      uicontrol(PM,'style','RadioButton',...
          'String','Yes',...
          'Position',[300 140 60 022].*WS,...
          'ToolTipString','Do Jacobian intensity modulation to compensate for compression/stretching',...
          'CallBack','FieldMap(''RadioButton'');',...
          'Tag','Jacobian',...
          'UserData',1);
      uicontrol(PM,'style','RadioButton',...
          'String','No',...
          'Position',[370 140 60 022].*WS,...
          'ToolTipString','Don''t do Jacobian intensity modulation to compensate for compression/stretching',...
          'CallBack','FieldMap(''RadioButton'');',...
          'Tag','Jacobian',...
          'UserData',2);
      
      uicontrol(PM,'Style','Edit',...
          'Position',[300 110 70 024].*WS,...
          'ToolTipString','Give total time for readout of EPI echo train ms',...
          'CallBack','FieldMap(''ReadOutTime'');',...
          'Tag','ReadTime');

    uicontrol(PM,'String','Load EPI image',...
          'Position',[30 70 120 022].*WS,...
          'ToolTipString','Load sample modulus EPI Analyze image',...
          'CallBack','FieldMap(''LoadEpi_Gui'');',...
          'Tag','LoadEpi');
      uicontrol(PM,'String','Match VDM to EPI',...
          'Position',[165 70 120 022].*WS,...
          'ToolTipString','Match vdm to EPI (but only if you want to...)',...
          'CallBack','FieldMap(''MatchVDM_gui'');',...
          'Tag','MatchVDM');
     uicontrol(PM,'String','Write unwarped',...
          'Position',[300 70 120 022].*WS,...
          'ToolTipString','Write unwarped EPI Analyze image',...
          'CallBack','FieldMap(''WriteUnwarped_Gui'');',...
          'Tag','WriteUnwarped');
     uicontrol(PM,'String','Load structural',...
          'Position',[30 30 120 022].*WS,...
          'ToolTipString','Load structural image (but only if you want to...)',...
          'CallBack','FieldMap(''LoadStructural_Gui'');',...
          'Tag','LoadStructural');
    uicontrol(PM,'String','Match structural',...
          'Position',[164 30 120 022].*WS,...
          'ToolTipString','Match structural image to EPI (but only if you want to...)',...
          'CallBack','FieldMap(''MatchStructural_Gui'');',...
          'Tag','MatchStructural');
    uicontrol(PM,'String','Help',...
          'Position',[320 30 80 022].*WS,...
          'ToolTipString','Help',...
          'CallBack','FieldMap(''Help_Gui'');',...
          'ForegroundColor','g',...
          'Tag','Help');
     uicontrol(PM,'String','Quit',...
          'Position',[420 30 60 022].*WS,...
          'Enable','On',...
          'ToolTipString','Exit toolbox',...
          'CallBack','FieldMap(''Quit'');',...
          'ForegroundColor','r',...
          'Tag','Quit');

     uicontrol(PM,'style','RadioButton',...
           'String','Yes',...
           'Position',[123 330 48 022].*WS,...
           'ToolTipString','Mask brain using magnitude image before processing',...
           'CallBack','FieldMap(''MaskBrain'');',...
           'Tag','MaskBrain',...
           'UserData',1);
     uicontrol(PM,'style','RadioButton',...
           'String','No',...
           'Position',[173 330 46 022].*WS,...
           'ToolTipString','Do not mask brain using magnitude image before processing',...
           'CallBack','FieldMap(''MaskBrain'');',...
           'Tag','MaskBrain',...
           'UserData',2);


% Find list of default files and set up menu accordingly

      def_files = FieldMap('GetDefaultFiles');
      menu_items = FieldMap('DefFiles2MenuItems',def_files);
      if length(menu_items)==1
         uicontrol(PM,'String',menu_items{1},...
           'Enable','Off',...
           'ToolTipString','Site specific default files',...
           'Position',[400 480 60 22].*WS);
      else
         uicontrol(PM,'Style','PopUp',...
           'String',menu_items,...
           'Enable','On',...
           'ToolTipString','Site specific default files',...
           'Position',[390 480 90 22].*WS,...
           'callback','FieldMap(''DefaultFile_Gui'');',...
           'UserData',def_files);          
      end
      
%-Apply defaults to buttons

      FieldMap('SynchroniseButtons');

      %
      % Disable necessary buttons and parameters until phase-map is finished.
      %
      FieldMap('Reset_Gui');
      
%=======================================================================
%
% Make sure buttons reflect current parameter values
%
%=======================================================================

   case 'synchronisebuttons'

    % Set input format

      if ~isempty(IP.uflags.iformat)
         if strcmp(IP.uflags.iformat,'RI');
            h = findobj(get(PM,'Children'),'Tag','IFormat','UserData',1);
        FieldMap('RadioOn',h);
            IP.uflags.iformat = 'RI';
         else
            h = findobj(get(PM,'Children'),'Tag','IFormat','UserData',2);
        FieldMap('RadioOn',h);
            IP.uflags.iformat = 'PM';
         end
      else % Default to RI
         h = findobj(get(PM,'Children'),'Tag','IFormat','UserData',1);
     FieldMap('RadioOn',h);
         IP.uflags.iformat = 'RI';
      end

      % Set brain masking defaults

      if ~isempty(IP.maskbrain)
         if IP.maskbrain == 1
            h = findobj(get(PM,'Children'),'Tag','MaskBrain','UserData',1);
            FieldMap('RadioOn',h);
            IP.maskbrain = 1;
         else
            h = findobj(get(PM,'Children'),'Tag','MaskBrain','UserData',2);
            FieldMap('RadioOn',h);
            IP.maskbrain = 0;
         end
      else % Default to no
         h = findobj(get(PM,'Children'),'Tag','MaskBrain','UserData',2);
         FieldMap('RadioOn',h);
         IP.maskbrain = 0;
      end
     
      % Set echo Times

      if ~isempty(IP.et{1})
         h = findobj(get(PM,'Children'),'Tag','ShortTime');
         set(h,'String',sprintf('%2.2f',IP.et{1}));
      end

      if ~isempty(IP.et{2})
         h = findobj(get(PM,'Children'),'Tag','LongTime');
         set(h,'String',sprintf('%2.2f',IP.et{2}));
      end


      % Set EPI field map

      if ~isempty(IP.epifm)
         if IP.epifm
            h = findobj(get(PM,'Children'),'Tag','EPI','UserData',1);
        FieldMap('RadioOn',h);
            IP.epifm = 1;
         else
            h = findobj(get(PM,'Children'),'Tag','EPI','UserData',2);
        FieldMap('RadioOn',h);
            IP.epifm = 0;
         end
      else % Default to yes
         h = findobj(get(PM,'Children'),'Tag','EPI','UserData',1);
     FieldMap('RadioOn',h);
         IP.epifm = 1;
      end

      % Set apply jacobian

      if ~isempty(IP.ajm)
         if IP.ajm
            h = findobj(get(PM,'Children'),'Tag','Jacobian','UserData',1);
        FieldMap('RadioOn',h);
            IP.ajm = 1;
         else
            h = findobj(get(PM,'Children'),'Tag','Jacobian','UserData',2);
        FieldMap('RadioOn',h);
            IP.ajm = 0;
         end
      else % Default to no
         h = findobj(get(PM,'Children'),'Tag','Jacobian','UserData',2);
     FieldMap('RadioOn',h);
         IP.ajm = 0;
      end
            
      % Set blip direction

      if ~isempty(IP.blipdir)
         if IP.blipdir == 1
            h = findobj(get(PM,'Children'),'Tag','BlipDir','UserData',1);
        FieldMap('RadioOn',h);
            IP.blipdir = 1;
         elseif IP.blipdir == -1
            h = findobj(get(PM,'Children'),'Tag','BlipDir','UserData',2);
        FieldMap('RadioOn',h);
            IP.blipdir = -1;
         else
            error('Invalid phase-encode direction');
         end
      else % Default to -1 = negative
         h = findobj(get(PM,'Children'),'Tag','BlipDir','UserData',1);
     FieldMap('RadioOn',h);
         IP.blipdir = -1;
      end

      % Set total readout time

      if ~isempty(IP.tert)
         h = findobj(get(PM,'Children'),'Tag','ReadTime');
         set(h,'String',sprintf('%2.2f',IP.tert));
      end

%=======================================================================
%
% Input format was changed 
%
%=======================================================================

   case 'inputformat'

      %
      % Enforce radio behaviour.
      %
      index = get(gcbo,'UserData');
      if index==1 
         partner = 2; 
         IP.uflags.iformat = 'RI';
      else 
         partner = 1; 
         IP.uflags.iformat = 'PM';
      end

      h = findobj(get(PM,'Children'),'Tag','IFormat','UserData',partner);
      maxv = get(gcbo,'Max');
      value = get(gcbo,'Value');

      if value == maxv
         set(h,'Value',get(h,'Min'));
      else
         set(h,'Value',get(h,'Max'));
      end 
            
      FieldMap('IFormat_Gui');
            
%=======================================================================
%
% A new default file has been chosen from the popup menu - gui version.
%
%=======================================================================

   case 'defaultfile_gui'
      m_file = get(gcbo,'UserData');
      m_file = m_file(get(gcbo,'Value'));
      m_file = m_file{1}(1:end-2);
      IP = FieldMap('SetParams',m_file);
      FieldMap('IFormat_Gui');
      FieldMap('SynchroniseButtons');
      
%=======================================================================
%
% Load real or imaginary part of dual echo-time data - gui version.
%
%=======================================================================

   case 'loadfile_gui'

      FieldMap('Reset_Gui');
      % 
      % File selection using gui
      %
      index = get(gcbo,'UserData');
      FieldMap('LoadFile',index);
      ypos = [445 420 385 360];
      uicontrol(PM,'Style','Text',...
          'String',spm_file(IP.P{index}.fname,['short',num2str(50)]),...
          'Position',[220 ypos(index) 260 20].*WS,...
          'ForegroundColor','k', 'BackgroundColor',col,...
          'FontName',PF.times, 'FontSize',FS(8));
      FieldMap('Assert');

%=======================================================================
%
% Load phase and magnitude images - gui version.
%
%=======================================================================

   case 'loadfilepm_gui'

      FieldMap('Reset_Gui');
      % 
      % File selection using gui
      %
      index = get(gcbo,'UserData');
      FieldMap('LoadFilePM',index);
      ypos = [445 420 385 360];
      uicontrol(PM,'Style','Text',...
          'String',spm_file(IP.P{index}.fname,['short',num2str(50)]),...
          'Position',[220 ypos(index) 260 20].*WS,...
          'ForegroundColor','k', 'BackgroundColor',col,...
          'FontName',PF.times, 'FontSize',FS(8));
      FieldMap('Assert');

%=======================================================================
%
% Create unwrapped phase-map - gui version
%
%=======================================================================

   case 'createfieldmap_gui'

      %
      % Create fieldmap from complex images...
      %

      status=FieldMap('CreateFieldMap',IP);
      if ~isempty(status)

         FieldMap('DisplayImage',FieldMap('MakedP'),[.05 .75 .95 .2],1);

         % Toggle relevant buttons
         FieldMap('ToggleGUI','Off','CreateFieldMap');
         FieldMap('ToggleGUI','On',char('LoadFieldMap','WriteFieldMap'));
         %
         % Check that have correct parameters ready before allowing
         % an image to be loaded for unwarping and hence conversion of 
         % fieldmap to voxel shift map.
         %
         if (FieldMap('GatherVDMData'))
            FieldMap('FM2VDM',IP);
            FieldMap('ToggleGUI','On',char('LoadEpi','EPI',...
                                           'BlipDir','Jacobian','ReadTime'));
         end
      end

%=======================================================================
%
% Load already prepared fieldmap (in Hz).
%
%=======================================================================

   case 'loadfieldmap_gui'

      %
      % Reset all previously loaded field maps and images to unwarp
      %
      FieldMap('Reset_gui');
      %
      % Select a precalculated phase map
      %

      FieldMap('LoadFieldMap');

      uicontrol(PM,'Style','Text',...
          'String',spm_file(IP.pP.fname,['short',num2str(50)]),...
          'Position',[220 307 260 20].*WS,...
          'ForegroundColor','k', 'BackgroundColor',col,...
          'FontName',PF.times, 'FontSize',FS(8));

      FieldMap('DisplayImage',FieldMap('MakedP'),[.05 .75 .95 .2],1);       
      FieldMap('ToggleGUI','Off',char('CreateFieldMap','WriteFieldMap',...
                                         'MatchVDM','WriteUnwarped',...
                                         'LoadStructural', 'MatchStructural'));
      
      % Check that have correct parameters ready before allowing
      % an image to be loaded for unwarping and hence conversion of 
      % fieldmap to voxel shift map.
      %
      if (FieldMap('GatherVDMData'))
         FieldMap('FM2VDM',IP);
         FieldMap('ToggleGUI','On',char('LoadEpi','EPI',...
                                           'BlipDir','Jacobian','ReadTime'));
      end

%=======================================================================
%
% Write out field map in Hz - using GUI
%
%=======================================================================

   case 'writefieldmap_gui'
  
      FieldMap('Write',IP.P{1},IP.fm.fpm,'fpm_',64,'Fitted phase map in Hz');
    
%=======================================================================
%
% Load sample EPI image using gui, unwarp and display result
%
%=======================================================================

   case 'loadepi_gui'

      %
      % Select and display EPI image
      %

      FieldMap('LoadEPI');
      FieldMap('DisplayImage',IP.epiP,[.05 .5 .95 .2],2);

      % The fm2vdm parameters may have been changed so must calculate
      % the voxel shift map with current parameters.
     
      if (FieldMap('GatherVDMData'))
         FieldMap('FM2VDM',IP);
         FieldMap('UnwarpEPI',IP);
         FieldMap('DisplayImage',IP.uepiP,[.05 .25 .95 .2],3);

         % Once EPI has been unwarped can enable unwarp checks and structural stuff
         FieldMap('ToggleGUI','On',char('MatchVDM','WriteUnwarped','LoadStructural'));
      end

      %
      % Redisplay phasemap in case .mat file associated with 
      % EPI image different from images that filadmap was
      % estimated from.
      %
      FieldMap('DisplayImage',FieldMap('MakedP'),[.05 .75 .95 .2],1);
      
%=======================================================================
%
% Coregister fieldmap magnitude image to EPI to unwarp using GUI 
%
%=======================================================================

   case 'matchvdm_gui'

   if (FieldMap('GatherVDMData'))
      FieldMap('MatchVDM',IP);
      FieldMap('DisplayImage',FieldMap('MakedP'),[.05 .75 .95 .2],1)
      FieldMap('UnwarpEPI',IP);
      FieldMap('DisplayImage',IP.epiP,[.05 .5 .95 .2],2);
      FieldMap('DisplayImage',IP.uepiP,[.05 .25 .95 .2],3);
      % Once EPI has been unwarped can enable unwarp checks
      FieldMap('ToggleGUI','On',char('MatchVDM','WriteUnwarped','LoadStructural'));
   end

%=======================================================================
%
% Load structural image
%
%=======================================================================

   case 'loadstructural_gui'

      FieldMap('LoadStructural');
      FieldMap('DisplayImage',IP.nwarp,[.05 0.0 .95 .2],4);
      %
      % Redisplay the other images to allow for recalculation
      % of size of display.
      %
      FieldMap('DisplayImage',FieldMap('MakedP'),[.05 .75 .95 .2],1);
      FieldMap('DisplayImage',IP.epiP,[.05 .5 .95 .2],2);
      FieldMap('DisplayImage',IP.uepiP,[.05 .25 .95 .2],3);
      
      FieldMap('ToggleGUI','On','MatchStructural');
      
%=======================================================================
%
% Match structural image to unwarped EPI
%
%=======================================================================

   case 'matchstructural_gui'

      FieldMap('MatchStructural',IP);
      FieldMap('DisplayImage',IP.nwarp,[.05 0.0 .95 .2],4);

%=======================================================================
%
% Write out unwarped EPI
%
%=======================================================================

   case 'writeunwarped_gui'
       
      unwarp_info=sprintf('Unwarped EPI:echo time difference=%2.2fms, EPI readout time=%2.2fms, Jacobian=%d',IP.uflags.etd, IP.tert,IP.ajm);    
      IP.uepiP = FieldMap('Write',IP.epiP,IP.uepiP.dat,'u',IP.epiP.dt(1),unwarp_info);
      FieldMap('ToggleGUI','On', 'LoadStructural');
      FieldMap('ToggleGUI','Off', 'WriteUnwarped');

%=======================================================================
%
% Help using spm_help
%
%=======================================================================

   case 'help_gui'

      FieldMap('help');

%=======================================================================
%
% Quit Toolbox
%
%=======================================================================

   case 'quit'

     Fgraph = spm_figure('FindWin','Graphics');
     if ~isempty(Fgraph)
        if DGW
            delete(Fgraph);
            DGW = 0;
        else
            spm_figure('Clear','Graphics');
        end    
      end
      PM = spm_figure('FindWin','FieldMap');
      if ~isempty(PM)
         delete(PM);
         PM = [];
      end

%=======================================================================
%
% The following functions are GUI related functions that go on behind 
% the scenes.
%=======================================================================
%
%=======================================================================
%
% Reset gui inputs when a new field map is to be calculated
%
%=======================================================================

   case 'reset_gui'

      % Clear any precalculated fieldmap name string
      uicontrol(PM,'Style','Text','String','',...
            'Position',[230 307 260 20].*WS,...
            'ForegroundColor','k','BackgroundColor',col);

      % Clear  Fieldmap values
      uicontrol(PM,'Style','Text',...
          'Position',[350 280 50 020].*WS,...
          'HorizontalAlignment','left',...
          'String','');

      % Disable input routes to make sure 
      % a new fieldmap is calculated.
      %
      FieldMap('ToggleGUI','Off',char('CreateFieldMap','WriteFieldMap',...
                                 'LoadEpi','WriteUnwarped','MatchVDM',...
                                 'LoadStructural','MatchStructural',...
                                 'EPI','BlipDir',...
                                 'Jacobian','ReadTime'));
      %
      % A new input image implies all processed data is void.
      %
      IP.fm = [];
      IP.vdm = [];
      IP.jim = [];
      IP.pP = [];
      IP.epiP = [];
      IP.uepiP = [];
      IP.vdmP = [];
      IP.fmagP = [];
      IP.wfmagP = [];
      ID = cell(4,1);

      spm_orthviews('Reset');
      spm_figure('Clear','Graphics');

%=======================================================================
%
% Use gui for Real and Imaginary or Phase and Magnitude.
%
%=======================================================================

   case 'iformat_gui'

      FieldMap('Reset_Gui');

      % Load Real and Imaginary or Phase and Magnitude Buttons
      if strcmp(IP.uflags.iformat,'PM')
         h=findobj('Tag','LoadRI');
         set(h,'Visible','Off');
         uicontrol(PM,'String','Load Phase',...
             'Position',[125 450 90 022].*WS,...
         'ToolTipString','Load Analyze image containing first phase image',...
         'CallBack','FieldMap(''LoadFilePM_Gui'');',...
             'UserData',1,...
             'Tag','LoadPM');
         uicontrol(PM,'String','Load Mag.',...
             'Position',[125 425 90 022].*WS,...
         'ToolTipString','Load Analyze image containing first magnitude image',...
         'CallBack','FieldMap(''LoadFilePM_Gui'');',...
             'UserData',2,...
             'Tag','LoadPM');
         uicontrol(PM,'String','Load Phase',...
             'Position',[125 390 90 022].*WS,...
         'ToolTipString','Load Analyze image containing second phase image',...
         'CallBack','FieldMap(''LoadFilePM_Gui'');',...
             'UserData',3,...
             'Tag','LoadPM');
         uicontrol(PM,'String','Load Mag.',...
             'Position',[125 365 90 022].*WS,...
         'ToolTipString','Load Analyze image containing second magnitudeimage',...
         'CallBack','FieldMap(''LoadFilePM_Gui'');',...
             'UserData',4,...
             'Tag','LoadPM');
      else

         % make the objects we don't want invisible
         h=findobj('Tag','LoadPM');
         set(h,'Visible','Off');

         uicontrol(PM,'String','Load Real',...
             'Position',[125 450 90 022].*WS,...
         'ToolTipString','Load Analyze image containing real part of short echo-time image',...
         'CallBack','FieldMap(''LoadFile_Gui'');',...
             'UserData',1,...
             'Tag','LoadRI');
         uicontrol(PM,'String','Load Imag',...
             'Position',[125 425 90 022].*WS,...
         'ToolTipString','Load Analyze image containing imaginary part of short echo-time image',...
         'CallBack','FieldMap(''LoadFile_Gui'');',...
             'UserData',2,...
             'Tag','LoadRI');
         uicontrol(PM,'String','Load Real',...
             'Position',[125 390 90 022].*WS,...
         'ToolTipString','Load Analyze image containing real part of long echo-time image',...
         'CallBack','FieldMap(''LoadFile_Gui'');',...
             'UserData',3,...
             'Tag','LoadRI');
         uicontrol(PM,'String','Load Imag',...
             'Position',[125 365 90 022].*WS,...
         'ToolTipString','Load Analyze image containing imaginary part of long echo-time image',...
         'CallBack','FieldMap(''LoadFile_Gui'');',...
             'UserData',4,...
             'Tag','LoadRI');
      end

%=======================================================================
%
% Brain masking or not
%
%=======================================================================

   case 'maskbrain'

      FieldMap('Reset_Gui');
      %
      % Enforce radio behaviour.
      %
      index = get(gcbo,'UserData');
      if index==1 
         partner=2; 
         IP.maskbrain=1;
      else 
         partner=1; 
         IP.maskbrain=0;
      end
      tag = get(gcbo,'Tag');
      value = get(gcbo,'Value');
      maxv = get(gcbo,'Max');
      h = findobj(get(PM,'Children'),'Tag',tag,'UserData',partner);
      if value == maxv
         set(h,'Value',get(h,'Min'));
      else
         set(h,'Value',get(h,'Max'));
      end 

      FieldMap('Assert');

%=======================================================================
%
% Echo time was changed or entered
%
%=======================================================================

   case 'echotime'

      FieldMap('Assert');

%=======================================================================
%
% Check if inputs are correct for calculating new phase map
%
%=======================================================================

   case 'assert'

      if ~isempty(IP.pP)  % If we're using precalculated fieldmap.
         go = 0;
      else
         go = 1;
         for i=1:2
            if isempty(IP.P{i}), go = 0; end
         end
         for i=3:4
            if (isempty(IP.P{i}) && strcmp(IP.uflags.iformat,'RI')), go = 0; end
         end  
         h = findobj(get(PM,'Children'),'Tag','ShortTime');
         IP.et{1} = str2num(get(h,'String'));
         h = findobj(get(PM,'Children'),'Tag','LongTime');
         IP.et{2} = str2num(get(h,'String'));
         if isempty(IP.et{1}) || isempty(IP.et{2}) || IP.et{2} < IP.et{1}
            go = 0;
         end
      end
      if go
         FieldMap('ToggleGui','On','CreateFieldMap');
      else % Switch off fieldmap creation
         FieldMap('ToggleGui','Off','CreateFieldMap');
      end   

%=======================================================================
%
% Enable or disable specified buttons or parameters in GUI
%
%=======================================================================

   case 'togglegui'

      h = get(PM,'Children');
      tags=varargin{3};
      for n=1:size(varargin{3},1)
         set(findobj(h,'Tag',deblank(tags(n,:))),'Enable',varargin{2});
      end

%=======================================================================
%
% A radiobutton was pressed.
%
%=======================================================================

   case 'radiobutton'

      %
      % Enforce radio behaviour.
      %
      index = get(gcbo,'UserData');
      if index==1 
         partner=2; 
      else 
         partner=1; 
      end
      tag = get(gcbo,'Tag');
      value = get(gcbo,'Value');
      maxv = get(gcbo,'Max');
      h = findobj(get(PM,'Children'),'Tag',tag,'UserData',partner);
      if value == maxv
         set(h,'Value',get(h,'Min'));
      else
         set(h,'Value',get(h,'Max'));
      end 

      % Update Hz to voxel displacement map if the input parameters are ok
      if (FieldMap('GatherVDMData'))

         FieldMap('FM2VDM',IP);
         %
         % Update unwarped EPI if one is loaded
         %
         if ~isempty(IP.epiP) 
            IP.uepiP = FieldMap('UnwarpEPI',IP);
            FieldMap('DisplayImage',IP.uepiP,[.05 .25 .95 .2],3);
            FieldMap('ToggleGUI','On',char('WriteUnwarped')); 
         end
      end

%=======================================================================
%
% Enforce radio-button behaviour
%
%=======================================================================

   case 'radioon'
      my_gcbo = varargin{2};
      index = get(my_gcbo,'UserData');
      if index==1 
         partner=2; 
      else 
         partner=1; 
      end
      set(my_gcbo,'Value',get(my_gcbo,'Max'));
      h = findobj(get(PM,'Children'),'Tag',get(my_gcbo,'Tag'),'UserData',partner);
      set(h,'Value',get(h,'Min'));      
      
%=======================================================================
%
% Total readout-time was changed or entered
%
%=======================================================================

   case 'readouttime'

      %
      % Update unwarped EPI
      %
      % This means the transformation phase-map to voxel displacement-map
      % has changed, which means the inversion will have to be redone,
      % which means (in the case of flash-based field map) coregistration
      % between field-map and sample EPI image should be redone. Phew!
      %
                                       
      if (FieldMap('GatherVDMData')) 
         FieldMap('FM2VDM',IP);
         % Update unwarped EPI if one is loaded
         if ~isempty(IP.epiP) 
            IP.uepiP = FieldMap('UnwarpEPI',IP);
            FieldMap('DisplayImage',IP.uepiP,[.05 .25 .95 .2],3);
            FieldMap('ToggleGUI','On',char('WriteUnwarped',...
                                          'LoadStructural'));
         end
         FieldMap('ToggleGUI','On',char('LoadEpi','EPI','BlipDir',...
                                           'Jacobian','ReadTime')); 
      else
         % If readtime is missing switch everything off...
         FieldMap('ToggleGUI','Off',char('LoadEpi','EPI',...
                                            'BlipDir','Jacobian',...
                                      'WriteUnwarped','LoadStructural',...
                                          'MatchStructural', 'MatchVDM'));
      end

%=======================================================================
%
% Sample UIControls pertaining to information needed for transformation
% phase-map -> voxel displacement-map
%
%=======================================================================

   case 'gathervdmdata'

      h = findobj(get(PM,'Children'),'Tag','EPI','UserData',1);
      v = get(h,{'Value','Max'});
      if v{1} == 1 
         IP.epifm = 1; 
      else
         IP.epifm = 0; 
      end
      h = findobj(get(PM,'Children'),'Tag','BlipDir','UserData',1);
      v = get(h,{'Value','Max'});
      if v{1} == 1 %% CHLOE: Changed 
         IP.blipdir = 1; 
      else
         IP.blipdir = -1; 
      end
      h = findobj(get(PM,'Children'),'Tag','Jacobian','UserData',1);
      v = get(h,{'Value','Max'});
      if v{1} == 1 %% CHLOE: Changed 
         IP.ajm = 1; 
      else
         IP.ajm = 0; 
      end
      h = findobj(get(PM,'Children'),'Tag','ReadTime');
      IP.tert = str2num(get(h,'String'));

      if isempty(IP.tert) || isempty(IP.fm) || isempty(IP.fm.fpm)       
         varargout{1} = 0;
         FieldMap('ToggleGui','On','ReadTime');
      else
         varargout{1} = 1;
      end
      
%=======================================================================
%
% Draw transversal and sagittal image.
%
%=======================================================================

   case 'redrawimages'
      global curpos;

      for i=1:length(ID)
         if ~isempty(ID{i}) && ~isempty(st.vols{i})
            set(st.vols{i}.ax{2}.ax,'Visible','Off'); % Disable event delivery in Coronal
            set(st.vols{i}.ax{2}.d,'Visible','Off');   % Make Coronal invisible
            set(st.vols{i}.ax{1}.ax,'Position',ID{i}.tra_pos); 
            set(st.vols{i}.ax{1}.ax,'ButtonDownFcn','FieldMap(''Orthviews'');');
            set(get(st.vols{i}.ax{1}.ax,'YLabel'),'String',ID{i}.label);      
            set(st.vols{i}.ax{3}.ax,'Position',ID{i}.sag_pos); 
            set(st.vols{i}.ax{3}.ax,'ButtonDownFcn','FieldMap(''Orthviews'');');
         end
      end
     
%=======================================================================
%
% Callback for orthviews
%
%=======================================================================

   case 'orthviews' 
   
   if strcmp(get(gcf,'SelectionType'),'normal')
      spm_orthviews('Reposition');
      %spm_orthviews('set_pos2cm');
   elseif strcmp(get(gcf,'SelectionType'),'extend')
      spm_orthviews('Reposition');
      %spm_orthviews('set_pos2cm');
      spm_orthviews('context_menu','ts',1);
   end;
   curpos = spm_orthviews('pos',1); 
   set(st.in, 'String',sprintf('%3.3f',spm_sample_vol(st.vols{1},curpos(1),curpos(2),curpos(3),st.hld)));

%=======================================================================
%
% Display image.
%
%=======================================================================

   case 'displayimage'

      Fgraph = spm_figure('FindWin','Graphics');
      
%       if isempty(Fgraph)
%          st.fig=spm_figure('Create','Graphics','Graphics');
%          DGW = 1;
%       end

    % Only open Graphics window if one has been found
      if isempty(Fgraph)
          return;
      end
          
      if ~isempty(ID{varargin{4}})
         spm_orthviews('Delete',ID{varargin{4}}.h);
         ID{varargin{4}} = [];
      end

      h = spm_orthviews('Image',varargin{2},[.01 .01 .01 .01]);

      % Set bounding box to allow display of all images
      spm_orthviews('MaxBB'); 
      
      % Put everything in space of original EPI image.
      if varargin{4} == 2
         %spm_orthviews('Space',h); % This was deleted for some reason
      end
      
      %
      % Get best possible positioning and scaling for
      % transversal and sagittal views.
      %
      tra_pos = get(st.vols{varargin{4}}.ax{1}.ax,'Position');
      sag_pos = get(st.vols{varargin{4}}.ax{3}.ax,'Position');
      field_pos = varargin{3};
      x_scale = field_pos(3) / (tra_pos(3) + sag_pos(3));
      height = max([tra_pos(4) sag_pos(4)]);
      y_scale = field_pos(4) / height;
      if x_scale > y_scale % Height limiting
         scale = y_scale;
         dx = (field_pos(3) - scale*(tra_pos(3) + sag_pos(3))) / 3;
         dy = 0;
      else % Width limiting
         scale = x_scale;
         dx = 0;
         dy = (field_pos(4) - scale*height) / 2;
      end
      tra_pos = [field_pos(1)+dx field_pos(2)+dy scale*tra_pos([3 4])];
      sag_pos = [field_pos(1)+tra_pos(3)+2*dx field_pos(2)+dy scale*sag_pos([3 4])];
      label = {'Fieldmap in Hz',...
                       'Warped EPI',...
                       'Unwarped EPI',...
                       'Structural'};
      ID{varargin{4}} = struct('tra_pos',    tra_pos,...
                               'sag_pos',    sag_pos,...
                               'h',          h,...
                               'label',      label{varargin{4}});
      FieldMap('RedrawImages');
    
      st.in = uicontrol(PM,'Style','Text',...
          'Position',[340 280 50 020].*WS,...
          'HorizontalAlignment','left',...
          'String','');
    
%=======================================================================
%=======================================================================
%
% The functions below are called by the various gui buttons but are
% not dependent on the gui to work. These functions can therefore also
% be called from a script bypassing the gui and can return updated 
% variables.
%
%=======================================================================
%=======================================================================
%
%=======================================================================
%
% Create and initialise parameters for FieldMap toolbox
%
%=======================================================================

   case 'initialise'
      
      %
      % Initialise parameters
      %
      ID = cell(4,1);
      IP.P = cell(1,4);
      IP.pP = [];
      IP.fmagP = [];
      IP.wfmagP = [];
      IP.epiP = [];
      IP.uepiP = [];
      IP.nwarp = [];
      IP.fm = [];
      IP.vdm = [];
      IP.et = cell(1,2);
      IP.epifm = [];
      IP.blipdir = [];
      IP.ajm = [];
      IP.tert = [];
      IP.vdmP = []; %% Check that this should be there %%
      IP.maskbrain = [];
      IP.uflags = struct('iformat','','method','','fwhm',[],'bmask',[],'pad',[],'etd',[],'ws',[]);
      IP.mflags = struct('template',[],'fwhm',[],'nerode',[],'ndilate',[],'thresh',[],'reg',[],'graphics',0);

      % Initially set brain mask to be empty 
      IP.uflags.bmask = [];

      % Set parameter values according to defaults
      FieldMap('SetParams');
      
      varargout{1}= IP;

%=======================================================================
%
% Sets parameters according to the defaults file that is being passed
%
%=======================================================================

   case 'setparams'
       if nargin == 1
           % "Default" default file
           pm_defaults;
       else
           % Scanner or sequence specific default file
           %eval(varargin{2});
           spm('Run',varargin{2});
           pm_def = spm('GetGlobal','pm_def');
       end

      % Define parameters for fieldmap creation
      IP.et{1} = pm_def.SHORT_ECHO_TIME;
      IP.et{2} = pm_def.LONG_ECHO_TIME;
      IP.maskbrain = pm_def.MASKBRAIN;

      % Set parameters for unwrapping
      IP.uflags.iformat = pm_def.INPUT_DATA_FORMAT;
      IP.uflags.method = pm_def.UNWRAPPING_METHOD;
      IP.uflags.fwhm = pm_def.FWHM;
      IP.uflags.pad = pm_def.PAD;
      IP.uflags.ws = pm_def.WS;
      IP.uflags.etd = pm_def.LONG_ECHO_TIME - pm_def.SHORT_ECHO_TIME;     

      % Set parameters for brain masking
      IP.mflags.template = pm_def.MFLAGS.TEMPLATE;
      IP.mflags.fwhm = pm_def.MFLAGS.FWHM;
      IP.mflags.nerode = pm_def.MFLAGS.NERODE;
      IP.mflags.ndilate = pm_def.MFLAGS.NDILATE;
      IP.mflags.thresh = pm_def.MFLAGS.THRESH;
      IP.mflags.reg = pm_def.MFLAGS.REG;
      IP.mflags.graphics = pm_def.MFLAGS.GRAPHICS;

      % Set parameters for unwarping 
      IP.ajm = pm_def.DO_JACOBIAN_MODULATION;
      IP.blipdir = pm_def.K_SPACE_TRAVERSAL_BLIP_DIR;
      IP.tert = pm_def.TOTAL_EPI_READOUT_TIME;
      IP.epifm = pm_def.EPI_BASED_FIELDMAPS;

      varargout{1}= IP;

%=======================================================================
%
% Get a list of all the default files that are present
%
%=======================================================================

   case 'getdefaultfiles'
      fmdir = fullfile(spm('Dir'),'toolbox','FieldMap');
      defnames = spm_select('FPList',fmdir,'^pm_defaults.*\.m');
      defnames = char(defnames,...
          spm_select('FPList',fullfile(fmdir,'FIL'),'^pm_defaults.*\.m'));
      varargout{1} = cellstr(defnames);

%=======================================================================
%
% Get list of menuitems from list of default files
%
%=======================================================================

   case 'deffiles2menuitems'
      for i=1:length(varargin{2})
         menit{i} = spm_file(varargin{2}{i},'basename');
         menit{i} = menit{i}(13:end); % "pm_defaults_<menit>.m"
         if isempty(menit{i}), menit{i} = 'Default'; end
      end
      varargout{1} = menit;
      
%=======================================================================
%
% Scale a phase map so that max = pi and min =-pi radians.
%
%=======================================================================

   case 'scale'

   F=varargin{2};
   V=spm_vol(F);
   vol=spm_read_vols(V);
   mn=min(vol(:));
   mx=max(vol(:));
   vol=-pi+(vol-mn)*2*pi/(mx-mn);
   V.dt(1)=4;   
   varargout{1} = FieldMap('Write',V,vol,'sc',V.dt(1),V.descrip);

%=======================================================================
%
% Load real or imaginary part of dual echo-time data - NO gui.
%
%=======================================================================

   case 'loadfile'

      index=varargin{2};
      prompt_string = {'Select short echo-time real',...
                       'Select short echo-time imaginary',...
                       'Select long echo-time real',...
                       'Select long echo-time imaginary'};
      %IP.P{index} = spm_vol(spm_get(1,'*.img',prompt_string{index}));
      % SPM
      IP.P{index} = spm_vol(spm_select(1,'image',prompt_string{index}));

      if isfield(IP,'pP') && ~isempty(IP.pP)
     IP.pP = [];
      end
      varargout{1} = IP.P{index};

%=======================================================================
%
% Load phase or magnitude part of dual (possibly) echo-time data - NO gui.
%
%=======================================================================

   case 'loadfilepm'

      index=varargin{2};
      prompt_string = {'Select phase image',...
                       'Select magnitude image',...
                       'Select second phase image or press done',...
                       'Select magnitude image or press done'};
      if index<3
         %IP.P{index} = spm_vol(spm_get(1,'*.img',prompt_string{index}));
         % SPM5 Update
         IP.P{index} = spm_vol(spm_select(1,'image',prompt_string{index}));
         if index==1
            doscl=spm_input('Scale this to radians?',1,'b','Yes|No',[1,0]);
            Finter=spm_figure('FindWin','Interactive');
            close(Finter);
            if doscl==1
               tmp=FieldMap('Scale',IP.P{index}.fname);
               IP.P{index} = spm_vol(tmp.fname);
            end
          end
      else
         %IP.P{index} = spm_vol(spm_get([0 1],'*.img',prompt_string{index}));
         % SPM5 Update
         IP.P{index} = spm_vol(spm_select([0 1],'image',prompt_string{index}));
         if index==3 && ~isempty(IP.P{index})
            doscl=spm_input('Scale this to radians?',1,'b','Yes|No',[1,0]);
            Finter=spm_figure('FindWin','Interactive');
            close(Finter);
            if doscl==1
               tmp=FieldMap('Scale',IP.P{index}.fname);
               IP.P{index} = spm_vol(tmp.fname);
            end
          end
      end
      if isfield(IP,'pP') && ~isempty(IP.pP)
     IP.pP = [];
      end
      varargout{1} = IP.P{index};

%=======================================================================
%
% Load already prepared fieldmap (in Hz) - no gui
%
%=======================================================================

   case 'loadfieldmap'

      %
      % Select field map
      %
      % 24/03/04 - Chloe change below to spm_get(1,'fpm_*.img')... 
      % IP.pP = spm_vol(spm_get(1,'fpm_*.img','Select field map'));
      % SPM5 Update
      IP.pP = spm_vol(spm_select(1,'^fpm.*\.img$','Select field map'));

      IP.fm.fpm = spm_read_vols(IP.pP);
      IP.fm.jac = pm_diff(IP.fm.fpm,2);
      if isfield(IP,'P') && ~isempty(IP.P{1})
     IP.P = cell(1,4);
      end
      varargout{1} = IP.fm;
      varargout{2} = IP.pP;

%=======================================================================
%
% Create unwrapped phase-map - NO gui
%
%=======================================================================

   case 'createfieldmap'

      IP=varargin{2};

      % First check that images are in same space. 

      if size([IP.P{1} IP.P{2} IP.P{3} IP.P{4}],2)==4
         ip_dim=cat(1,[IP.P{1}.dim' IP.P{2}.dim' IP.P{3}.dim' IP.P{4}.dim']');
         %ip_mat=cat(2,[IP.P{1}.mat(:) IP.P{2}.mat(:) IP.P{3}.mat(:) IP.P{4}.mat(:)]');
         ip_mat=cat(2,single([IP.P{1}.mat(:) IP.P{2}.mat(:) IP.P{3}.mat(:) IP.P{4}.mat(:)]'));
      else
         ip_dim=cat(1,[IP.P{1}.dim' IP.P{2}.dim']');
         %ip_mat=cat(2,[IP.P{1}.mat(:) IP.P{2}.mat(:)]');
         ip_mat=cat(2,single([IP.P{1}.mat(:) IP.P{2}.mat(:)]'));
      end

      if any(any(diff(ip_dim,1,1),1)&[1,1,1])
        errordlg({'Images don''t all have same dimensions'});
        drawnow;
        varargout{1}=[];
      elseif any(any(abs(diff(ip_mat,1,1))>1e-4))
        errordlg({'Images don''t all have same orientation & voxel size'});
        drawnow;
        varargout{1}=[];
      else
         % Update flags for unwarping (in case TEs have been adjusted
         IP.uflags.etd = IP.et{2}-IP.et{1};      

         % Clear any brain mask 
         IP.uflags.bmask = [];

         % SPM5 Update
         % If flag selected to mask brain and the field map is not based on EPI
         if IP.maskbrain==1 
            IP.fmagP = FieldMap('Magnitude',IP);
            IP.uflags.bmask = pm_brain_mask(IP.fmagP,IP.mflags);
            varargout{2} = IP.fmagP;
         end

         IP.fm = pm_make_fieldmap([IP.P{1} IP.P{2} IP.P{3} IP.P{4}],IP.uflags);
         varargout{1} = IP.fm;
      end

%=======================================================================
%
% Convert field map to voxel displacement map and 
% do necessary inversions of displacement fields.
%
%=======================================================================

   case 'fm2vdm'

      IP=varargin{2};
      %
      % If we already have memory mapped image pointer to voxel
      % displacement map let's reuse it (so not to void possible
      % realignment). If field-map is non-EPI based the previous
      % realignment (based on different parameters) will be non-
      % optimal and we should advice user to redo it.
      %
      % If no pointer already exist we'll make one.
      %
      if isfield(IP,'vdmP') && ~isempty(IP.vdmP)
     msgbox({'Changing this parameter means that if previously',...
                 'you matched VDM to EPI, this result may no longer',...
                 'be optimal. In this case we recommend you redo the',...
                 'Match VDM to EPI.'},...
         'Coregistration notice','warn','modal');
      else
         if ~isempty(IP.pP)
            IP.vdmP = struct('dim',    IP.pP.dim(1:3),...
                             'dt',[4 spm_platform('bigend')],...
                             'mat',    IP.pP.mat);
            IP.vdmP.fname = spm_file(IP.pP.fname,'prefix','vdm5_');
         else
            IP.vdmP = struct('dim',    IP.P{1}.dim(1:3),...
                             'dt',[4 spm_platform('bigend')],...
                             'mat',    IP.P{1}.mat);
            IP.vdmP.fname = spm_file(IP.P{1}.fname,'prefix','vdm5_');
         end
      end

      % Scale field map and jacobian by total EPI readout time
      IP.vdm.vdm = IP.blipdir*IP.tert*1e-3*IP.fm.fpm;
      IP.vdm.jac = IP.blipdir*IP.tert*1e-3*IP.fm.jac;

      % Chloe added this: 26/02/05
      % Put fieldmap parameters in descrip field of vdm

      vdm_info=sprintf('Voxel Displacement Map:echo time difference=%2.2fms, EPI readout time=%2.2fms',IP.uflags.etd, IP.tert);    
      if IP.epifm ==1
         spm_progress_bar('Init',3,'Inverting displacement map','');
         spm_progress_bar('Set',1);
         % Invert voxel shift map and multiply by -1...
         IP.vdm.ivdm = pm_invert_phasemap(-1*IP.vdm.vdm);
         IP.vdm.ijac = pm_diff(IP.vdm.ivdm,2); 
         spm_progress_bar('Set',2);
         spm_progress_bar('Clear');
         FieldMap('Write',IP.vdmP,IP.vdm.ivdm,'',IP.vdmP.dt(1),vdm_info);
      else
         FieldMap('Write',IP.vdmP,IP.vdm.vdm,'',IP.vdmP.dt(1),vdm_info);
      end
      varargout{1} = IP.vdm;
      varargout{2} = IP.vdmP;

%=======================================================================
%
% Write out image
%
%=======================================================================

   case 'write'
  
      V    = varargin{2};
      vol  = varargin{3};
      name = varargin{4};

      % Write out image
      img = struct(...
          'fname', spm_file(V.fname,'prefix',name),...
          'dim',    V.dim(1:3),...
          'dt',     [varargin{5} spm_platform('bigend')],...
          'mat',    V.mat,...
          'descrip',varargin{6});
      
      img = spm_write_vol(img,vol);
      varargout{1} = img;

%=======================================================================
%
% Create fieldmap (Hz) struct for use when displaying image.
%
%=======================================================================

   case 'makedp'

   if isfield(IP,'vdmP') && ~isempty(IP.vdmP);
      dP = struct('dim',  IP.vdmP.dim,...
                  'dt',[64 spm_platform('bigend')],...
                  'pinfo',   [1 0]',...
                  'mat',     IP.vdmP.mat,...
                  'dat',     reshape(IP.fm.fpm,IP.vdmP.dim),...
                  'fname',   'display_image');
   else
      if isfield(IP,'P') && ~isempty(IP.P{1})
         dP = struct('dim',     IP.P{1}.dim,...
                     'dt',[64 spm_platform('bigend')],...
                     'pinfo',   [1 0]',...
                     'mat',     IP.P{1}.mat,...
                     'dat',     reshape(IP.fm.fpm,IP.P{1}.dim),...
                     'fname',   'display_image');
      elseif isfield(IP,'pP') && ~isempty(IP.pP)
         dP = struct('dim',     IP.pP.dim,...
                     'dt',[64 spm_platform('bigend')],...
                     'pinfo',   [1 0]',...
                     'mat',     IP.pP.mat,...
                     'dat',     reshape(IP.fm.fpm,IP.pP.dim),...
                     'fname',   'display_image');
      end
   end
   varargout{1} = dP;

%=======================================================================
%
% Load sample EPI image - NO gui
%
%=======================================================================

   case 'loadepi'

   %
   % Select EPI
   %
   %IP.epiP = spm_vol(spm_get(1,'*.img','Select sample EPI image'));
   % SPM5 Update
   IP.epiP = spm_vol(spm_select(1,'image','Select sample EPI image'));

   varargout{1} = IP.epiP;

%=======================================================================
%
% Create unwarped epi - NO gui
%
%=======================================================================

   case 'unwarpepi'
     
      %
      % Update unwarped EPI
      %
      IP=varargin{2};
      IP.uepiP = struct('fname',   'Image in memory',...
                        'dim',     IP.epiP.dim,...
                        'dt',[64 spm_platform('bigend')],...
                        'pinfo',   IP.epiP.pinfo(1:2),...
                        'mat',     IP.epiP.mat);

      % Need to sample EPI and voxel shift map in space of EPI...
      [x,y,z] = ndgrid(1:IP.epiP.dim(1),1:IP.epiP.dim(2),1:IP.epiP.dim(3));
      xyz = [x(:) y(:) z(:)];

      % Space of EPI is IP.epiP{1}.mat and space of 
      % voxel shift map is IP.vdmP{1}.mat 
      tM = inv(IP.epiP.mat\IP.vdmP.mat);

      x2 = tM(1,1)*x + tM(1,2)*y + tM(1,3)*z + tM(1,4);
      y2 = tM(2,1)*x + tM(2,2)*y + tM(2,3)*z + tM(2,4);
      z2 = tM(3,1)*x + tM(3,2)*y + tM(3,3)*z + tM(3,4);
      xyz2 = [x2(:) y2(:) z2(:)];

      %
      % Make mask since it is only meaningful to calculate undistorted
      % image in areas where we have information about distortions.
      %
      msk = reshape(double(xyz2(:,1)>=1 & xyz2(:,1)<=IP.vdmP.dim(1) &...
                   xyz2(:,2)>=1 & xyz2(:,2)<=IP.vdmP.dim(2) &...
                   xyz2(:,3)>=1 & xyz2(:,3)<=IP.vdmP.dim(3)),IP.epiP.dim(1:3));
              
      % Read in voxel displacement map in correct space
      tvdm = reshape(spm_sample_vol(spm_vol(IP.vdmP.fname),xyz2(:,1),...
                      xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));

      % Voxel shift map must be added to the y-coordinates. 
      uepi = reshape(spm_sample_vol(IP.epiP,xyz(:,1),...
                      xyz(:,2)+tvdm(:),xyz(:,3),1),IP.epiP.dim(1:3));% TEMP CHANGE
      
      % Sample Jacobian in correct space and apply if required 
      if IP.ajm==1
         if IP.epifm==1 % If EPI, use inverted jacobian

            IP.jim = reshape(spm_sample_vol(IP.vdm.ijac,xyz2(:,1),...
                      xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));
         else
            IP.jim = reshape(spm_sample_vol(IP.vdm.jac,xyz2(:,1),...
                      xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));
         end  
         uepi = uepi.*(1+IP.jim);
      end

      IP.uepiP.dat=uepi.*msk;
      varargout{1}=IP.uepiP;
      
%=======================================================================
%
% Create unwarped epi - NO gui ***FOR XY***
%
%=======================================================================

   case 'unwarpepixy'
      %
      % Update unwarped EPI
      %
      IP=varargin{2};
      IP.uepiP = struct('fname',   'Image in memory',...
                        'dim',     IP.epiP.dim,...
                        'dt',[64 spm_platform('bigend')],...
                        'pinfo',   IP.epiP.pinfo(1:2),...
                        'mat',     IP.epiP.mat);

      % Need to sample EPI and voxel shift map in space of EPI...
      [x,y,z] = ndgrid(1:IP.epiP.dim(1),1:IP.epiP.dim(2),1:IP.epiP.dim(3));
      xyz = [x(:) y(:) z(:)];

      % Space of EPI is IP.epiP{1}.mat and space of 
      % voxel shift map is IP.vdmP{1}.mat 
      tM = inv(IP.epiP.mat\IP.vdmP.mat);

      x2 = tM(1,1)*x + tM(1,2)*y + tM(1,3)*z + tM(1,4);
      y2 = tM(2,1)*x + tM(2,2)*y + tM(2,3)*z + tM(2,4);
      z2 = tM(3,1)*x + tM(3,2)*y + tM(3,3)*z + tM(3,4);
      xyz2 = [x2(:) y2(:) z2(:)];

      %
      % Make mask since it is only meaningful to calculate undistorted
      % image in areas where we have information about distortions.
      %
      msk = reshape(double(xyz2(:,1)>=1 & xyz2(:,1)<=IP.vdmP.dim(1) &...
                           xyz2(:,2)>=1 & xyz2(:,2)<=IP.vdmP.dim(2) &...
                           xyz2(:,3)>=1 & xyz2(:,3)<=IP.vdmP.dim(3)),IP.epiP.dim(1:3));
              
      % Read in voxel displacement map in correct space
      tvdm = reshape(spm_sample_vol(spm_vol(IP.vdmP.fname),xyz2(:,1),...
                      xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));

      % Voxel shift map must be added to the x-coordinates. 
      uepi = reshape(spm_sample_vol(IP.epiP,xyz(:,1)+tvdm(:),...
                      xyz(:,2),xyz(:,3),1),IP.epiP.dim(1:3));% TEMP CHANGE
      
      % Sample Jacobian in correct space and apply if required 
      if IP.ajm==1
         if IP.epifm==1 % If EPI, use inverted jacobian

            IP.jim = reshape(spm_sample_vol(IP.vdm.ijac,xyz2(:,1),...
                      xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));
         else
            IP.jim = reshape(spm_sample_vol(IP.vdm.jac,xyz2(:,1),...
                      xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));
         end  
         uepi = uepi.*(1+IP.jim);
      end

      IP.uepiP.dat=uepi.*msk;
      varargout{1}=IP.uepiP;

%=======================================================================
%
% Coregister fieldmap magnitude image to EPI to unwarp 
%
%=======================================================================

   case 'matchvdm'
   
      IP=varargin{2};

      % 
      % Need a fieldmap magnitude image
      %

      if isempty(IP.pP) && ~isempty(IP.P{1})

         IP.fmagP = struct(...
            'fname', spm_file(IP.P{1}.fname,'prefix','mag_'),...
            'dim',   IP.P{1}.dim,...
            'dt',    IP.P{1}.dt,...
            'pinfo', IP.P{1}.pinfo,...
            'mat',   IP.P{1}.mat);
      
         % If using real and imaginary data, calculate using sqrt(i1.^2 + i2.^2).
         % If using phase and magnitude, use magnitude image.
         if strcmp(IP.uflags.iformat,'RI') 
            IP.fmagP = spm_imcalc(spm_vol([IP.P{1}.fname;IP.P{2}.fname]),IP.fmagP,'sqrt(i1.^2 + i2.^2)');
         else
            IP.fmagP = IP.P{2};
         end
      elseif ~isempty(IP.pP) && ~isempty(IP.fmagP)
          msg=sprintf('Using %s for matching\n',IP.fmagP.fname);
          disp(msg);
      else
         IP.fmagP = spm_vol(spm_select(1,'image','Select field map magnitude image'));
      end

      % Now we have field map magnitude image, we want to coregister it to the 
      % EPI to be unwarped. 
      % If using an EPI field map:
      % 1) Coregister magnitude image to EPI.
      % 2) Apply resulting transformation matrix to voxel shift map
      % If using a non-EPI field map:
      % 1) Forward warp magnitude image
      % 2) Coregister warped magnitude image to EPI.
      % 3) Apply resulting transformation matrix to voxel shift map
    
      if IP.epifm==1
         [mi,M] = FieldMap('Coregister',IP.epiP,IP.fmagP);
         MM = IP.fmagP.mat;
      else
         % Need to sample magnitude image in space of EPI to be unwarped...
         [x,y,z] = ndgrid(1:IP.epiP.dim(1),1:IP.epiP.dim(2),1:IP.epiP.dim(3));
          xyz = [x(:) y(:) z(:)];

         % Space of EPI is IP.epiP{1}.mat and space of fmagP is IP.fmagP.mat 
         tM = inv(IP.epiP.mat\IP.fmagP.mat);
         x2 = tM(1,1)*x + tM(1,2)*y + tM(1,3)*z + tM(1,4);
         y2 = tM(2,1)*x + tM(2,2)*y + tM(2,3)*z + tM(2,4);
         z2 = tM(3,1)*x + tM(3,2)*y + tM(3,3)*z + tM(3,4);
         xyz2 = [x2(:) y2(:) z2(:)];
         wfmag = reshape(spm_sample_vol(IP.fmagP,xyz2(:,1),...
                      xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));     

         % Need to sample voxel shift map in space of EPI to be unwarped
         tvdm = reshape(spm_sample_vol(IP.vdm.vdm,xyz2(:,1),...
                      xyz2(:,2),xyz2(:,3),0),IP.epiP.dim(1:3));            
                      
         % Now apply warps to resampled forward warped magnitude image...
         wfmag = reshape(spm_sample_vol(wfmag,xyz(:,1),xyz(:,2)-tvdm(:),...
                       xyz(:,3),1),IP.epiP.dim(1:3));
              
         % Write out forward warped magnitude image
         IP.wfmagP = struct('dim',  IP.epiP.dim,...
                         'dt',[64 spm_platform('bigend')],...
                         'pinfo',   IP.epiP.pinfo,...
                         'mat',     IP.epiP.mat);
         IP.wfmagP = FieldMap('Write',IP.epiP,wfmag,'wfmag_',4,'Voxel shift map');

         % Now coregister warped magnitude field map to EPI
         [mi,M] = FieldMap('Coregister',IP.epiP,IP.wfmagP);

         % Update the .mat file of the forward warped mag image 
         spm_get_space(deblank(IP.wfmagP.fname),M*IP.wfmagP.mat);

         % Get the original space of the fmap magnitude      
         MM = IP.fmagP.mat;
      end

      % Update .mat file for voxel displacement map
      IP.vdmP.mat=M*MM;       
      spm_get_space(deblank(IP.vdmP.fname),M*MM); 

      varargout{1} = IP.vdmP; 
      
%=======================================================================
%
% Coregister fieldmap magnitude image to EPI to do unwarpxy
%
%=======================================================================

   case 'matchvdmxy'
   
      IP=varargin{2};

      % 
      % Need a fieldmap magnitude image
      %

      if isempty(IP.pP) && ~isempty(IP.P{1})

         IP.fmagP=struct(...
             'fname', spm_file(IP.P{1}.fname,'prefix','mag_'),...
             'dim',   IP.P{1}.dim,...
             'dt',    IP.P{1}.dt,...
             'pinfo', IP.P{1}.pinfo,...
             'mat',   IP.P{1}.mat);
      
         % If using real and imaginary data, calculate using sqrt(i1.^2 + i2.^2).
         % If using phase and magnitude, use magnitude image.
         if strcmp(IP.uflags.iformat,'RI') 
            IP.fmagP = spm_imcalc(spm_vol([IP.P{1}.fname;IP.P{2}.fname]),IP.fmagP,'sqrt(i1.^2 + i2.^2)');
         else
            IP.fmagP = IP.P{2};
         end
      elseif ~isempty(IP.pP) && ~isempty(IP.fmagP)
          msg=sprintf('Using %s for matching\n',IP.fmagP.fname);
          disp(msg);
      else
         %IP.fmagP = spm_vol(spm_get(1,'*.img','Select field map magnitude image'));
         % SPM5 Update
         IP.fmagP = spm_vol(spm_select(1,'image','Select field map magnitude image'));
      end

      % Now we have field map magnitude image, we want to coregister it to the 
      % EPI to be unwarped. 
      % If using an EPI field map:
      % 1) Coregister magnitude image to EPI.
      % 2) Apply resulting transformation matrix to voxel shift map
      % If using a non-EPI field map:
      % 1) Forward warp magnitude image
      % 2) Coregister warped magnitude image to EPI.
      % 3) Apply resulting transformation matrix to voxel shift map
    
      if IP.epifm==1
         [mi,M] = FieldMap('Coregister',IP.epiP,IP.fmagP);
         MM = IP.fmagP.mat;
      else
         % Need to sample magnitude image in space of EPI to be unwarped...
         [x,y,z] = ndgrid(1:IP.epiP.dim(1),1:IP.epiP.dim(2),1:IP.epiP.dim(3));
          xyz = [x(:) y(:) z(:)];

         % Space of EPI is IP.epiP{1}.mat and space of fmagP is IP.fmagP.mat 
         tM = inv(IP.epiP.mat\IP.fmagP.mat);
         x2 = tM(1,1)*x + tM(1,2)*y + tM(1,3)*z + tM(1,4);
         y2 = tM(2,1)*x + tM(2,2)*y + tM(2,3)*z + tM(2,4);
         z2 = tM(3,1)*x + tM(3,2)*y + tM(3,3)*z + tM(3,4);
         xyz2 = [x2(:) y2(:) z2(:)];
         wfmag = reshape(spm_sample_vol(IP.fmagP,xyz2(:,1),...
                      xyz2(:,2),xyz2(:,3),1),IP.epiP.dim(1:3));     

         % Need to sample voxel shift map in space of EPI to be unwarped
         tvdm = reshape(spm_sample_vol(IP.vdm.vdm,xyz2(:,1),...
                      xyz2(:,2),xyz2(:,3),0),IP.epiP.dim(1:3));            
                      
         % Now apply warps to resampled forward warped magnitude image...
         wfmag = reshape(spm_sample_vol(wfmag,xyz(:,1)-tvdm(:),xyz(:,2),...
                       xyz(:,3),1),IP.epiP.dim(1:3));
              
         % Write out forward warped magnitude image
         IP.wfmagP = struct('dim',  IP.epiP.dim,...
                         'dt',[64 spm_platform('bigend')],...
                         'pinfo',   IP.epiP.pinfo,...
                         'mat',     IP.epiP.mat);
         IP.wfmagP = FieldMap('Write',IP.epiP,wfmag,'wfmag_',4,'Voxel shift map');

         % Now coregister warped magnitude field map to EPI
         [mi,M] = FieldMap('Coregister',IP.epiP,IP.wfmagP);

         % Update the .mat file of the forward warped mag image 
         spm_get_space(deblank(IP.wfmagP.fname),M*IP.wfmagP.mat);

         % Get the original space of the fmap magnitude      
         MM = IP.fmagP.mat;
      end

      % Update .mat file for voxel displacement map
      IP.vdmP.mat=M*MM;       
      spm_get_space(deblank(IP.vdmP.fname),M*MM); 

      varargout{1} = IP.vdmP; 
     
%=======================================================================
%
% Invert voxel displacement map
%
%=======================================================================

   case 'invert'

   vdm = pm_invert_phasemap(IP.vdm.vdm);
   varargout{1} = vdm; 

%=======================================================================
%
% Get fieldmap magnitude image
%
%=======================================================================

   case 'magnitude'

   IP=varargin{2};

   if isempty(IP.fmagP)
      % 
      % Get fieldmap magnitude image
      %
      if isempty(IP.pP) && ~isempty(IP.P{1})

         IP.fmagP=struct(...
             'fname', spm_file(IP.P{1}.fname,'prefix','mag_'),...
             'dim',   IP.P{1}.dim,...
             'dt',    IP.P{1}.dt,...                         
             'pinfo', IP.P{1}.pinfo,...
             'mat',   IP.P{1}.mat);
      
         % If using real and imaginary data, calculate using sqrt(i1.^2 + i2.^2).
         % If using phase and magnitude, use magnitude image.
         if IP.uflags.iformat=='RI' 

               IP.fmagP = spm_imcalc(spm_vol([IP.P{1}.fname;IP.P{2}.fname]),IP.fmagP,'sqrt(i1.^2 + i2.^2)');
         else
               IP.fmagP = IP.P{2};
         end

     % else
     %    IP.fmagP = spm_vol(spm_get(1,'*.img','Select field map magnitude image'));

      else
       do_f = ~isfield(IP, 'fmagP');
       if ~do_f, do_f = isempty(IP.fmagP); end
       if do_f
         IP.fmagP = spm_vol(spm_select(1,'image','Select field map magnitude image'));
       end

      end
   end 
   varargout{1} = IP.fmagP;
      
%=======================================================================
%
% Load structural image
%
%=======================================================================

case 'loadstructural'
      %IP.nwarp = spm_vol(spm_get(1,'*.img','Select structural image'));
      % SPM5 Update
      IP.nwarp = spm_vol(spm_select(1,'image','Select structural image'));

      varargout{1} = IP.nwarp;

%=======================================================================
%
% Coregister structural image to unwarped EPI
%
%=======================================================================

   case 'matchstructural'
   
     IP = varargin{2};
     [mi,M] = FieldMap('Coregister',IP.uepiP,IP.nwarp);
     MM = IP.nwarp.mat;
     IP.nwarp.mat=M*MM;
     spm_get_space(deblank(IP.nwarp.fname),IP.nwarp.mat);
     varargout{1}=mi;

%=======================================================================
%
% Coregister two images - NO gui
%
%=======================================================================

   case 'coregister'
     
      %
      % Coregisters image VF to image VG (use VF and VG like in SPM code)
      %
      VG = varargin{2};
      VF = varargin{3};

      % Define flags for coregistration...
      IP.cflags.cost_fun = spm_get_defaults('coreg.estimate.cost_fun');
      IP.cflags.sep      = spm_get_defaults('coreg.estimate.sep');
      IP.cflags.tol      = spm_get_defaults('coreg.estimate.tol');
      IP.cflags.fwhm     = spm_get_defaults('coreg.estimate.fwhm'); 
      IP.cflags.graphics = 0;

      % Voxel sizes (mm)
      vxg   = sqrt(sum(VG.mat(1:3,1:3).^2));
      vxf   = sqrt(sum(VF.mat(1:3,1:3).^2));

      % Smoothnesses
      fwhmg = sqrt(max([1 1 1]*IP.cflags.sep(end)^2 - vxg.^2, [0 0 0]))./vxg;
      fwhmf = sqrt(max([1 1 1]*IP.cflags.sep(end)^2 - vxf.^2, [0 0 0]))./vxf;

   
      % Need to load data smoothed in unit8 format (as in spm_coreg_ui)
%      if isfield(VG,'dat')
%         VG=rmfield(VG,'dat');
%      end
      if ~isfield(VG,'uint8')
         VG.uint8 = loaduint8(VG);
     VG       = smooth_uint8(VG,fwhmg); % Note side effects
      end

 %     if isfield(VF,'dat')
 %        VF=rmfield(VF,'dat');
 %     end
      if ~isfield(VF,'uint8')
         VF.uint8 = loaduint8(VF);
         VF       = smooth_uint8(VF,fwhmf); % Note side effects
      end

      x=spm_coreg(VG,VF,IP.cflags);
      M  = inv(spm_matrix(x));
      MM = spm_get_space(deblank(VF.fname));

      varargout{1}=spm_coreg(x,VG,VF,2,IP.cflags.cost_fun,IP.cflags.fwhm);
      varargout{2}=M;

%=======================================================================
%
% Use spm_help to display help for FieldMap toolbox
%
%=======================================================================

   case 'help'
      spm_help('FieldMap.man');
   return
end

%=======================================================================
%
% Smoothing functions required for spm_coreg
%
%=======================================================================

function V = smooth_uint8(V,fwhm)

% Convolve the volume in memory (fwhm in voxels).
lim = ceil(2*fwhm);
s  = fwhm/sqrt(8*log(2));
x  = [-lim(1):lim(1)]; x = smoothing_kernel(fwhm(1),x); x  = x/sum(x);
y  = [-lim(2):lim(2)]; y = smoothing_kernel(fwhm(2),y); y  = y/sum(y);
z  = [-lim(3):lim(3)]; z = smoothing_kernel(fwhm(3),z); z  = z/sum(z);
i  = (length(x) - 1)/2;
j  = (length(y) - 1)/2;
k  = (length(z) - 1)/2;
spm_conv_vol(V.uint8,V.uint8,x,y,z,-[i j k]);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function krn = smoothing_kernel(fwhm,x)

% Variance from FWHM
s = (fwhm/sqrt(8*log(2)))^2+eps;

% The simple way to do it. Not good for small FWHM
% krn = (1/sqrt(2*pi*s))*exp(-(x.^2)/(2*s));

% For smoothing images, one should really convolve a Gaussian
% with a sinc function.  For smoothing histograms, the
% kernel should be a Gaussian convolved with the histogram
% basis function used. This function returns a Gaussian
% convolved with a triangular (1st degree B-spline) basis
% function.

% Gaussian convolved with 0th degree B-spline
% int(exp(-((x+t))^2/(2*s))/sqrt(2*pi*s),t= -0.5..0.5)
% w1  = 1/sqrt(2*s);
% krn = 0.5*(erf(w1*(x+0.5))-erf(w1*(x-0.5)));

% Gaussian convolved with 1st degree B-spline
%  int((1-t)*exp(-((x+t))^2/(2*s))/sqrt(2*pi*s),t= 0..1)
% +int((t+1)*exp(-((x+t))^2/(2*s))/sqrt(2*pi*s),t=-1..0)
w1  =  0.5*sqrt(2/s);
w2  = -0.5/s;
w3  = sqrt(s/2/pi);
krn = 0.5*(erf(w1*(x+1)).*(x+1) + erf(w1*(x-1)).*(x-1) - 2*erf(w1*x   ).* x)...
      +w3*(exp(w2*(x+1).^2)     + exp(w2*(x-1).^2)     - 2*exp(w2*x.^2));

krn(krn<0) = 0;
return;

%=======================================================================
%
% Load data function required for spm_coreg
%
%=======================================================================

function udat = loaduint8(V)
% Load data from file indicated by V into an array of unsigned bytes.
if size(V.pinfo,2)==1 && V.pinfo(1) == 2
    mx = 255*V.pinfo(1) + V.pinfo(2);
    mn = V.pinfo(2);
else
    spm_progress_bar('Init',V.dim(3),...
        ['Computing max/min of ' spm_file(V.fname,'filename')],...
        'Planes complete');
    mx = -Inf; mn =  Inf;
    for p=1:V.dim(3)
        img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
        mx  = max([max(img(:))+paccuracy(V,p) mx]);
        mn  = min([min(img(:)) mn]);
        spm_progress_bar('Set',p);
    end
end
spm_progress_bar('Init',V.dim(3),...
    ['Loading ' spm_file(V.fname,'filename')],...
    'Planes loaded');

udat = uint8(0);
udat(V.dim(1),V.dim(2),V.dim(3))=0;
rand('state',100);
for p=1:V.dim(3)
    img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
    acc = paccuracy(V,p);
    if acc==0
        udat(:,:,p) = uint8(round((img-mn)*(255/(mx-mn))));
    else
        % Add random numbers before rounding to reduce aliasing artifact
        r = rand(size(img))*acc;
        udat(:,:,p) = uint8(round((img+r-mn)*(255/(mx-mn))));
    end
    spm_progress_bar('Set',p);
end
spm_progress_bar('Clear');
return;

function acc = paccuracy(V,p)
if ~spm_type(V.dt(1),'intt')
    acc = 0;
else
    if size(V.pinfo,2)==1
        acc = abs(V.pinfo(1,1));
    else
        acc = abs(V.pinfo(1,p));
    end
end
%_______________________________________________________________________
