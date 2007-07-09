function [status] = hastoolbox(toolbox, add_to_path);

% HASTOOLBOX tests whether an external toolbox is installed. Optionally
% it will try to determine the path to the toolbox and install it
% automatically.
% 
% Use as
%   [status] = hastoolbox(toolbox, add_to_path);

% Copyright (C) 2005-2006, Robert Oostenveld
%
% $Log: hastoolbox.m,v $
% Revision 1.10  2007/05/06 09:10:07  roboos
% added spm5
%
% Revision 1.9  2007/02/26 13:41:07  roboos
% made small change to fastica detection (suggested by Sameer)
%
% Revision 1.8  2007/02/13 17:22:27  roboos
% added MRI from eeg.sf.net
%
% Revision 1.7  2007/02/13 14:01:26  roboos
% added brainstorm
%
% Revision 1.6  2007/02/12 19:43:23  roboos
% added fastica, optim
%
% Revision 1.5  2007/01/17 17:05:34  roboos
% added matlab signal processing toolbox
%
% Revision 1.4  2007/01/04 12:25:19  roboos
% added SON2
%
% Revision 1.3  2007/01/03 17:01:15  roboos
% added 4d-version toolbox
%
% Revision 1.2  2006/06/07 10:48:02  roboos
% changed the "see xxx" string
%
% Revision 1.1  2006/06/07 09:28:41  roboos
% renamed fieldtrip/private/checktoolbox into misc/hastoolbox
%
% Revision 1.8  2006/06/06 14:18:22  roboos
% added neuroshare, eeprobe, yokogawa
%
% Revision 1.7  2006/05/31 08:56:24  roboos
% implemented detection of toolbox in users ~/matlab/toolboxname
%
% Revision 1.6  2006/05/23 10:20:46  roboos
% added beowulf and mentat toolboxes
%
% Revision 1.5  2006/04/26 11:37:22  roboos
% added besa toolbox
%
% Revision 1.4  2006/02/07 20:01:39  roboos
% aded biosig and meg-pd (neuromag)
%
% Revision 1.3  2006/01/17 14:05:54  roboos
% added GLNA64 for mentat000
%
% Revision 1.2  2006/01/06 11:39:23  roboos
% added copyrigth and cvs logging, changed some comments
%

% this points the user to the website where he/she can download the toolbox
url = {
  'AFNI'       'see http://afni.nimh.nih.gov'
  'DSS'        'see http://www.cis.hut.fi/projects/dss'
  'EEGLAB'     'see http://www.sccn.ucsd.edu/eeglab'
  'NWAY'       'see http://www.models.kvl.dk/source/nwaytoolbox'
  'SPM2'       'see http://www.fil.ion.ucl.ac.uk/spm'
  'SPM5'       'see http://www.fil.ion.ucl.ac.uk/spm'
  'MEG-PD'     'see http://www.kolumbus.fi/kuutela/programs/meg-pd'
  'MEG-CALC'   'this is a commercial toolbox from Neuromag, see http://www.neuromag.com'
  'BIOSIG'     'see http://biosig.sourceforge.net'
  'EEG'        'see http://eeg.sourceforge.net'
  'EEGSF'      'see http://eeg.sourceforge.net'  % alternative name
  'MRI'        'see http://eeg.sourceforge.net'  % alternative name
  'NEUROSHARE' 'see http://www.neuroshare.org'
  'BESA'       'see http://www.megis.de, or contact Karsten Hoechstetter'
  'EEPROBE'    'see http://www.ant-neuro.com, or contact Maarten van der Velde'
  'YOKOGAWA'   'see http://www.yokogawa.co.jp, or contact Nobuhiko Takahashi'
  'BEOWULF'    'see http://oase.uci.ru.nl/~roberto, or contact Robert Oostenveld'
  'MENTAT'     'see http://oase.uci.ru.nl/~roberto, or contact Robert Oostenveld'
  'SON2'       'see http://www.kcl.ac.uk/depsta/biomedical/cfnr/lidierth.html, or contact Malcolm Lidierth' 
  '4D-VERSION' 'contact Christian Wienbruch'
  'SIGNAL'     'see http://www.mathworks.com/products/signal'
  'OPTIM'      'see http://www.mathworks.com/products/optim'
  'FASTICA'    'see http://www.cis.hut.fi/projects/ica/fastica'
  'BRAINSTORM' 'see http://neuroimage.ucs.edu/brainstorm'
};

if nargin<2
  % default is not to add the path automatically
  add_to_path = 0;
end

% determine whether the toolbox is installed
toolbox = upper(toolbox);
switch toolbox
  case 'AFNI'
    status = (exist('BrikLoad') && exist('BrikInfo'));
  case 'DSS'
    status = exist('dss', 'file') && exist('dss_create_state', 'file');
  case 'EEGLAB'
    status = exist('runica', 'file');
  case 'NWAY'
    status = exist('parafac', 'file');
  case 'SPM2'
    status = exist('spm_vol') && exist('spm_write_vol') && exist('spm_normalise');
  case 'SPM5'
    status = exist('spm_vol') && exist('spm_write_vol') && exist('spm_normalise') && exist('spm_vol_nifti');
  case 'MEG-PD'
    status = (exist('rawdata') && exist('channames'));
  case 'MEG-CALC'
    status = (exist('megmodel') && exist('megfield') && exist('megtrans'));
  case 'BIOSIG'
    status = (exist('sopen') && exist('sread'));
  case 'EEG'
    status = (exist('ctf_read_res4') && exist('ctf_read_meg4'));
  case 'EEGSF'  % alternative name
    status = (exist('ctf_read_res4') && exist('ctf_read_meg4'));
  case 'MRI'    % other functions in the mri section
    status = (exist('avw_hdr_read') && exist('avw_img_read'));
  case 'NEUROSHARE'
    status  = (exist('ns_OpenFile') && exist('ns_SetLibrary') && exist('ns_GetAnalogData'));
  case 'BESA'
    status = (exist('readBESAtfc') && exist('readBESAswf'));
  case 'EEPROBE'
    status  = (exist('read_eep_avr') && exist('read_eep_cnt'));
  case 'YOKOGAWA'
    status  = (exist('GetMeg160ChannelInfoM') && exist('GetMeg160ContinuousRawDataM'));
  case 'BEOWULF'
    status = (exist('evalwulf') && exist('evalwulf') && exist('evalwulf'));
  case 'MENTAT'
    status  = (exist('pcompile') && exist('pfor') && exist('peval'));
  case 'SON2'
    status  = (exist('SONFileHeader') && exist('SONChanList') && exist('SONGetChannel'));
  case '4D-VERSION'
    status  = (exist('read4d') && exist('read4dhdr'));
  case 'SIGNAL'
    status = exist('medfilt1');
  case 'OPTIM'
    status  = (exist('fmincon') && exist('fminunc'));
  case 'FASTICA'
    status  = exist('fastica', 'file');
  case 'BRAINSTORM'
    status  = exist('bem_xfer');
  otherwise
    warning(sprintf('cannot determine whether the %s toolbox is present', toolbox));
    status = 0;
end
% it should be a boolean value
status = (status~=0);

if ~status && add_to_path 
  % try to determine the path of the toolbox
  if (strcmp(computer, 'GLNX86') || strcmp(computer, 'GLNXA64')) && isdir('/home/common/matlab/')
    % for linux computers in the F.C. Donders Centre
    prefix = '/home/common/matlab/';
  elseif strcmp(computer, 'PCWIN') && isdir('h:\common\matlab\')
    % for windows computers in the F.C. Donders Centre
    prefix = 'h:\common\matlab\';
  elseif isdir([getenv('HOME') '/matlab/'])
    % use the matlab subdirectory in your homedirectory (works on unix and mac)
    prefix = [getenv('HOME') '/matlab/'];
  else
    prefix = [];
  end
  toolboxpath = [prefix lower(toolbox)];
  if exist(toolboxpath)
    % add the toolbox to the path automatically
    warning(sprintf('adding %s toolbox to your Matlab path', toolbox));
    addpath(toolboxpath);
    status = 1;
  else
    % the toolbox is not on the path and cannot be added
    sel = find(strcmp(url(:,1), toolbox));
    if ~isempty(sel)
      msg = sprintf('the %s toolbox is not installed, %s', toolbox, url{sel, 2});
    else
      msg = sprintf('the %s toolbox is not installed', toolbox);
    end
    error(msg);
  end
end

