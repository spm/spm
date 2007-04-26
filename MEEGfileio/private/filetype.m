function [ftype, detail] = filetype(filename, desired);

% FILETYPE determines the filetype of many EEG/MEG/MRI data files by
% looking at the name, extension and optionally (part of) its contents.
% It tries to determine the global type of file (which usually
% corresponds to the manufacturer, the recording system or to the
% software used to create the file) and the particular subtype (e.g.
% continuous, average).
%
% Use as
%   type = filetype(filename)
%   type = filetype(dirname)
% This gives you a descriptive string with the data type, and can be
% used in a switch-statement. The descriptive string that is returned
% usually is something like 'XXX_YYY'/ where XXX refers to the
% manufacturer and YYY to the type of the data.
%
% Alternatively, use as
%   flag = filetype(filename, type)
%   flag = filetype(dirname, type)
% This gives you a boolean flag (0 or 1) indicating whether the file
% is of the desired type, and can be used to check whether the
% user-supplied file is what your subsequent code expects.
%
% Alternatively, use as
%   flag = filetype(dirlist, type)
% where the dirlist contains a list of files contained within one
% directory. This gives you a boolean vector indicating for each file
% whether it is of the desired type.
%
% Most filetypes of the following manufacturers and/or software programs are recognized
%  - VSMMedtech/CTF
%  - Elektra/Neuromag
%  - Yokogawa
%  - 4D/BTi
%  - EDF
%  - Neuroscan
%  - Analyse
%  - EEProbe
%  - BrainVision
%  - BESA
%  - Curry
%  - ASA
%  - LORETA
%  - Analyze/SPM
%  - MINC
%  - AFNI
%  - Neuralynx
%  - Plexon
%
% See also READ_XXX_YYY where XXX=manufacturer and YYY=subtype

% Copyright (C) 2003-2006 Robert Oostenveld
%
% $Log: filetype.m,v $
% Revision 1.45  2007/03/21 17:21:30  roboos
% remove . and .. from the file listing in case of a directory as input
% removed neuralynx_nte, the correct extension is *.nts
% added header check to neuralynx_nts
% implemented subfunction most, c.f. any
% swiched from using any(...) to most(...) for determining content of dataset directory
% implemented plexon_ds for directory with nex files in it
% made some additional small changes
%
% Revision 1.44  2007/03/19 16:52:37  roboos
% added neuralynx_nte
%
% Revision 1.43  2007/01/09 09:29:25  roboos
% small change
%
% Revision 1.42  2007/01/04 08:12:00  roboos
% fixed bug for besa_avr, renamed an incorrect tag into plexon_plx
%
% Revision 1.41  2006/12/12 21:00:56  roboos
% fixed bug for neuralynx_cds
%
% Revision 1.40  2006/12/12 11:50:19  roboos
% added support for complete listing of directories, implemented detection of neuralynx_sdma using recursion
%
% Revision 1.39  2006/11/30 10:05:22  roboos
% reinserted support for Yokogawa sqd, as yokogawa_ave
%
% Revision 1.38  2006/09/18 21:51:55  roboos
% implemented support for fcdc_matbin, i.e. a dataset consisting of a matlab file with header and events and a seperate binary datafile
%
% Revision 1.37  2006/09/18 14:23:50  roboos
% added detection of 4D-BTi file formats
%
% Revision 1.36  2006/08/28 10:11:19  roboos
% moved subfunctions check_header and check_extension into seperate files and renamed them
% changed some | into ||
%
% Revision 1.35  2006/06/07 15:09:07  roboos
% improved the check for matlab, look at the file header
%
% Revision 1.34  2006/04/26 08:56:17  roboos
% added detection of neuralynx_sdma
%
% Revision 1.33  2006/04/24 12:18:36  roboos
% added split DMA log file
%
% Revision 1.32  2006/04/19 07:53:43  roboos
% added *.ima for DICOM
%
% Revision 1.31  2006/03/30 07:40:11  roboos
% added loreta, cleaned up documentation
%
% Revision 1.30  2006/03/23 16:42:50  roboos
% implemented support for neuralynx nse, nts, nvt
% implemented support for neuralynx_cds, i.e. a combined dataset with seperate directories for lfp, mua and spike channels
%
% Revision 1.29  2006/03/16 17:32:57  roboos
% added besa_swf
%
% Revision 1.28  2006/03/14 09:14:33  roboos
% added neuralynx_tsl, neuralynx_tsh and neuralynx_ttl
%
% Revision 1.27  2006/03/09 12:33:37  roboos
% added support for multiple files in check_extension subfunction
% improved detection of neuralynx_ds, also when Events.Nev is absent
%
% Revision 1.26  2006/02/24 12:50:33  roboos
% changed strcmp(lower(..)) into strcmpi(..)
%
% Revision 1.25  2005/12/16 14:02:20  roboos
% added (commented out) line for bham_bdf
%
% Revision 1.24  2005/12/02 08:38:09  roboos
% added biosemi_bdf EEG data format
%
% Revision 1.23  2005/11/30 11:37:07  roboos
% added support for ced_son
%
% Revision 1.22  2005/09/29 00:48:01  roboos
% dropped *.sqd and changed from yokogawa_sqd to yokogawa_ave (as suggested by Masahiro Shimogawara)
% added support for MBFYS tri and ama files
%
% Revision 1.21  2005/10/04 16:52:24  roboos
% added support for MPI Frankfurt *.dap files
% added check for directory in check_header
% moved directory listing in case of directory to the beginning of the function
% some small typo's and help changes
%
% Revision 1.20  2005/09/09 09:13:32  roboos
% added support for neuralynx_dma
%
% Revision 1.19  2005/09/06 12:43:41  roboos
% renamed plextor into plexon (incorrect company name)
%
% Revision 1.18  2005/09/02 16:28:00  roboos
% changed Plexon *.ddt from plexon_nex to plexon_ddt, since they are different formats and not compatible with each other
%
% Revision 1.17  2005/09/01 09:24:01  roboos
% added Yokogawa MEG formats
% added BESA beamfourmer source reconstruction (besa_src)
% prevent type-clash between brainvision_dat and besa_src
%
% Revision 1.16  2005/08/04 07:43:49  roboos
% added neuralynx log file
%
% Revision 1.15  2005/07/29 11:27:11  roboos
% added BESA tfc, mul and gen
%
% Revision 1.14  2005/05/18 06:49:52  roboos
% changed single & to double &&
%
% Revision 1.13  2005/05/12 07:22:40  roboos
% added Neuralynx filetypes and dataset directory
%
% Revision 1.12  2005/04/27 06:18:57  roboos
% added support for MEG42RS as a valid CTF res4 subformat
%
% Revision 1.11  2005/04/08 06:47:44  roboos
% added BESA paradigm file *.pdg
%
% Revision 1.10  2005/02/16 08:06:32  roboos
% changed behavior of magic number check: if file cannot be read, give warning (instead of error) and indicate false
% changed indentation
%
% Revision 1.9  2005/02/11 07:53:35  roboos
% iadded support for Curry, Polhemus, Plexon and  Besa *.eps
%
% Revision 1.8  2005/02/07 19:26:19  roboos
% added Ole Jensen's *.lay fileformat
%
% Revision 1.7  2004/10/28 09:42:19  roboos
% added AFNI brick and head (based on extension)
%
% Revision 1.6  2004/08/26 15:57:09  roboos
% changed file extension check into case insensitive
% added support for brainvision_marker
%
% Revision 1.5  2004/08/26 11:35:30  roboos
% added MINC based on mnc extension
%
% Revision 1.4  2004/03/29 15:14:30  roberto
% added support for *.dat as brainvision_dat (exported from BVA)
%
% Revision 1.3  2004/03/23 14:46:13  roberto
% fixed bug in eeprobe average, made besa average smarter
%
% Revision 1.2  2004/03/10 14:06:41  roberto
% added Neuromag filetypes, updated and improved help
%
% Revision 1.1  2003/12/05 10:59:43  roberto
% initial submission
%

if nargin<2
  desired = [];
end

if iscell(filename)
  % perform the test for each filename, return a boolean vector
  ftype = false(size(filename));
  for i=1:length(filename)
    if strcmp(filename{i}(end), '.')
      % do not recurse into this directory or the parent directory
      continue
    else
      ftype(i) = filetype(filename{i}, desired);
    end
  end
  return
end

% start with unknown values
ftype        = 'unknown';
manufacturer = 'unknown';
content      = 'unknown';

[p, f, x] = fileparts(filename);

if isdir(filename)
  % the directory listing is needed below
  ls = dir(filename);
  % remove the parent directory and the directory itself from the list
  ls = ls(~strcmp({ls.name}, '.'));
  ls = ls(~strcmp({ls.name}, '..'));
  for i=1:length(ls)
    % make sure that the directory listing includes the complete path
    ls(i).name = fullfile(filename, ls(i).name);
  end
end

% known CTF file types
if filetype_check_extension(filename, '.ds') && isdir(filename)
  ftype = 'ctf_ds';
  manufacturer = 'CTF';
  content = 'MEG dataset';
elseif filetype_check_extension(filename, '.res4') && (filetype_check_header(filename, 'MEG41RS') || filetype_check_header(filename, 'MEG42RS') || filetype_check_header(filename, 'MEG4RES'))
  ftype = 'ctf_res4';
  manufacturer = 'CTF';
  content = 'MEG/EEG header information';
elseif filetype_check_extension(filename, '.meg4') && filetype_check_header(filename, 'MEG41CP')
  ftype = 'ctf_meg4';
  manufacturer = 'CTF';
  content = 'MEG/EEG';
elseif filetype_check_extension(filename, '.mrk') && filetype_check_header(filename, 'PATH OF DATASET:')
  ftype = 'ctf_mrk';
  manufacturer = 'CTF';
  content = 'marker file';
elseif filetype_check_extension(filename, '.mri') && filetype_check_header(filename, 'CTF_MRI_FORMAT VER 2.2')
  ftype = 'ctf_mri';
  manufacturer = 'CTF';
  content = 'MRI';
elseif filetype_check_extension(filename, '.hdm')
  ftype = 'ctf_hdm';
  manufacturer = 'CTF';
  content = 'volume conduction model';
elseif filetype_check_extension(filename, '.hc')
  ftype = 'ctf_hc';
  manufacturer = 'CTF';
  content = 'headcoil locations';
elseif filetype_check_extension(filename, '.shape')
  ftype = 'ctf_shape';
  manufacturer = 'CTF';
  content = 'headshape points';
elseif filetype_check_extension(filename, '.shape_info')
  ftype = 'ctf_shapeinfo';
  manufacturer = 'CTF';
  content = 'headshape information';

  % known Neuromag file types
elseif filetype_check_extension(filename, '.fif')
  ftype = 'neuromag_fif';
  manufacturer = 'Neuromag';
  content = 'MEG header and data';
elseif filetype_check_extension(filename, '.bdip')
  ftype = 'neuromag_bdip';
  manufacturer = 'Neuromag';
  content = 'dipole model';

  % known Yokogawa file types
elseif filetype_check_extension(filename, '.ave') || filetype_check_extension(filename, '.sqd')
  ftype = 'yokogawa_ave';
  manufacturer = 'Yokogawa';
  content = 'averaged MEG data';
elseif filetype_check_extension(filename, '.con')
  ftype = 'yokogawa_con';
  manufacturer = 'Yokogawa';
  content = 'continuous MEG data';
elseif filetype_check_extension(filename, '.raw')
  ftype = 'yokogawa_raw';
  manufacturer = 'Yokogawa';
  content = 'evoked/trialbased MEG data';
elseif filetype_check_extension(filename, '.mri') && filetype_check_header(filename, '####')  % FIXME, not correct
  ftype = 'yokogawa_mri';
  manufacturer = 'Yokogawa';
  content = 'anatomical MRI';

elseif filetype_check_header(filename, 'E|lk') % I am not sure whether this always applies
  ftype = '4d_pdf';
  manufacturer = '4D/BTI';
  content = 'raw MEG data (processed data file)';
elseif exist([filename '.m4d'], 'file') && exist([filename '.xyz'], 'file') % these two ascii header files accompany the raw data
  ftype = '4d_pdf';
  manufacturer = '4D/BTI';
  content = 'raw MEG data (processed data file)';
elseif filetype_check_extension(filename, '.m4d') && exist([filename(1:(end-3)) 'xyz'], 'file') % these come in pairs
  ftype = '4d_m4d';
  manufacturer = '4D/BTI';
  content = 'MEG header information';
elseif filetype_check_extension(filename, '.xyz') && exist([filename(1:(end-3)) 'm4d'], 'file') % these come in pairs
  ftype = '4d_xyz';
  manufacturer = '4D/BTI';
  content = 'MEG sensor positions';

  % known EEProbe file types
elseif filetype_check_extension(filename, '.cnt') && filetype_check_header(filename, 'RIFF')
  ftype = 'eep_cnt';
  manufacturer = 'EEProbe';
  content = 'EEG';
elseif filetype_check_extension(filename, '.avr') && filetype_check_header(filename, char([38 0 16 0]))
  ftype = 'eep_avr';
  manufacturer = 'EEProbe';
  content = 'ERP';
elseif filetype_check_extension(filename, '.trg')
  ftype = 'eep_trg';
  manufacturer = 'EEProbe';
  content = 'trigger information';
elseif filetype_check_extension(filename, '.rej')
  ftype = 'eep_rej';
  manufacturer = 'EEProbe';
  content = 'rejection marks';

  % known ASA file types
elseif filetype_check_extension(filename, '.elc')
  ftype = 'asa_elc';
  manufacturer = 'ASA';
  content = 'electrode positions';
elseif filetype_check_extension(filename, '.vol')
  ftype = 'asa_vol';
  manufacturer = 'ASA';
  content = 'volume conduction model';
elseif filetype_check_extension(filename, '.bnd')
  ftype = 'asa_bnd';
  manufacturer = 'ASA';
  content = 'boundary element model details';
elseif filetype_check_extension(filename, '.msm')
  ftype = 'asa_msm';
  manufacturer = 'ASA';
  content = 'ERP';
elseif filetype_check_extension(filename, '.msr')
  ftype = 'asa_msr';
  manufacturer = 'ASA';
  content = 'ERP';
elseif filetype_check_extension(filename, '.dip')
  % FIXME, can also be CTF dipole file
  ftype = 'asa_dip';
  manufacturer = 'ASA';
elseif filetype_check_extension(filename, '.mri')
  % FIXME, can also be CTF mri file
  ftype = 'asa_mri';
  manufacturer = 'ASA';
  content = 'MRI image header';
elseif filetype_check_extension(filename, '.iso')
  ftype = 'asa_iso';
  manufacturer = 'ASA';
  content = 'MRI image data';

  % known Neuroscan file types
elseif filetype_check_extension(filename, '.avg') && filetype_check_header(filename, 'Version 3.0')
  ftype = 'ns_avg';
  manufacturer = 'Neuroscan';
  content = 'averaged EEG';
elseif filetype_check_extension(filename, '.cnt') && filetype_check_header(filename, 'Version 3.0')
  ftype = 'ns_cnt';
  manufacturer = 'Neuroscan';
  content = 'continuous EEG';
elseif filetype_check_extension(filename, '.eeg') && filetype_check_header(filename, 'Version 3.0')
  ftype = 'ns_eeg';
  manufacturer = 'Neuroscan';
  content = 'epoched EEG';

  % known Analyze & SPM file types
elseif filetype_check_extension(filename, '.hdr')
  ftype = 'analyze_hdr';
  manufacturer = 'Mayo Analyze';
  content = 'PET/MRI image header';
elseif filetype_check_extension(filename, '.img')
  ftype = 'analyze_img';
  manufacturer = 'Mayo Analyze';
  content = 'PET/MRI image data';
elseif filetype_check_extension(filename, '.mnc')
  ftype = 'minc';
  content = 'MRI image data';

  % known LORETA file types
elseif filetype_check_extension(filename, '.lorb')
  ftype = 'loreta_lorb';
  manufacturer = 'old LORETA';
  content = 'source reconstruction';
elseif filetype_check_extension(filename, '.slor')
  ftype = 'loreta_slor';
  manufacturer = 'sLORETA';
  content = 'source reconstruction';


  % known AFNI file types
elseif filetype_check_extension(filename, '.brik') || filetype_check_extension(filename, '.BRIK')
  ftype = 'afni_brik';
  content = 'MRI image data';
elseif filetype_check_extension(filename, '.head') || filetype_check_extension(filename, '.HEAD')
  ftype = 'afni_head';
  content = 'MRI header data';

  % known BrainVison file types
elseif filetype_check_extension(filename, '.vhdr')
  ftype = 'brainvision_vhdr';
  manufacturer = 'BrainProducts';
  content = 'EEG header';
elseif filetype_check_extension(filename, '.vmrk')
  ftype = 'brainvision_vmrk';
  manufacturer = 'BrainProducts';
  content = 'EEG markers';
elseif filetype_check_extension(filename, '.vabs')
  ftype = 'brainvision_vabs';
  manufacturer = 'BrainProducts';
  content = 'Brain Vison Analyzer macro';
elseif filetype_check_extension(filename, '.eeg')
  % FIXME, can also be Neuroscan epoched EEG data
  ftype = 'brainvision_eeg';
  manufacturer = 'BrainProducts';
  content = 'continuous EEG data';
elseif filetype_check_extension(filename, '.seg')
  ftype = 'brainvision_seg';
  manufacturer = 'BrainProducts';
  content = 'segmented EEG data';
elseif filetype_check_extension(filename, '.dat') && ~filetype_check_header(filename, 'BESA_SA_IMAGE')
  % WARNING this is a very general name, it could be exported BrainVision
  % data but also a BESA beamformer source reconstruction
  ftype = 'brainvision_dat';
  manufacturer = 'BrainProducts';
  content = 'exported EEG data';
elseif filetype_check_extension(filename, '.marker')
  ftype = 'brainvision_marker';
  manufacturer = 'BrainProducts';
  content = 'rejection markers';

  % known Polhemus file types
elseif filetype_check_extension(filename, '.pos')
  ftype = 'polhemus_pos';
  manufacturer = 'BrainProducts/CTF/Polhemus?'; % actually I don't know whose software it is
  content = 'electrode positions';

  % known Neuralynx file types
elseif filetype_check_extension(filename, '.nev')
  ftype = 'neuralynx_nev';
  manufacturer = 'Neuralynx';
  content = 'event information';
elseif filetype_check_extension(filename, '.ncs') && filetype_check_header(filename, '####')
  ftype = 'neuralynx_ncs';
  manufacturer = 'Neuralynx';
  content = 'continuous single channel recordings';
elseif filetype_check_extension(filename, '.nse') && filetype_check_header(filename, '####')
  ftype = 'neuralynx_nse';
  manufacturer = 'Neuralynx';
  content = 'spike waveforms';
elseif filetype_check_extension(filename, '.nts')  && filetype_check_header(filename, '####')
  ftype = 'neuralynx_nts';
  manufacturer = 'Neuralynx';
  content = 'timestamps only';
elseif filetype_check_extension(filename, '.nvt')
  ftype = 'neuralynx_nvt';
  manufacturer = 'Neuralynx';
  content = 'video tracker';
elseif filetype_check_extension(filename, '.ntt')
  ftype = 'neuralynx_ntt';
  manufacturer = 'Neuralynx';
  content = 'continuous tetrode recordings';
elseif strcmpi(f, 'logfile') && strcmpi(x, '.txt')  % case insensitive
  ftype = 'neuralynx_log';
  manufacturer = 'Neuralynx';
  content = 'log information in ASCII format';
elseif ~isempty(strfind(lower(f), 'dma')) && strcmpi(x, '.log')  % this is not a very strong detection
  ftype = 'neuralynx_dma';
  manufacturer = 'Neuralynx';
  content = 'raw aplifier data directly from DMA';
elseif isdir(filename) && ~isempty(strmatch('Events.Nev', {ls.name}))
  % a regular Neuralynx dataset directory contains an event file
  ftype = 'neuralynx_ds';
  manufacturer = 'Neuralynx';
  content = 'dataset';
elseif isdir(filename) && most(filetype_check_extension({ls.name}, '.ncs'))
  % a directory containing continuously sampled channels in Neuralynx format
  ftype = 'neuralynx_ds';
  manufacturer = 'Neuralynx';
  content = 'continuously sampled channels';
elseif isdir(filename) && most(filetype_check_extension({ls.name}, '.nse'))
  % a directory containing spike waveforms in Neuralynx format
  ftype = 'neuralynx_ds';
  manufacturer = 'Neuralynx';
  content = 'spike waveforms';
elseif isdir(filename) && most(filetype_check_extension({ls.name}, '.nte'))
  % a directory containing spike timestamps in Neuralynx format
  ftype = 'neuralynx_ds';
  manufacturer = 'Neuralynx';
  content = 'spike timestamps';

  % these are formally not Neuralynx file formats, but at the FCDC we use them together with Neuralynx
elseif isdir(filename) && any(filetype_check_extension({ls.name}, '.ttl')) && any(filetype_check_extension({ls.name}, '.crc'))
  % a directory containing the split channels from a DMA logfile
  ftype = 'neuralynx_sdma';
  manufacturer = 'F.C. Donders Centre';
  content = 'dataset';
elseif isdir(filename) && any(filetype({ls.name}, 'neuralynx_ds'))
  % a downsampled Neuralynx DMA file can be split into three seperate lfp/mua/spike directories
  % treat them as one combined dataset
  ftype = 'neuralynx_cds';
  manufacturer = 'F.C. Donders Centre';
  content = 'dataset containing seperate lfp/mua/spike directories';
elseif filetype_check_extension(filename, '.tsl') && filetype_check_header(filename, 'tsl')
  ftype = 'neuralynx_tsl';
  manufacturer = 'F.C. Donders Centre';
  content = 'timestamps from DMA log file';
elseif filetype_check_extension(filename, '.tsh') && filetype_check_header(filename, 'tsh')
  ftype = 'neuralynx_tsh';
  manufacturer = 'F.C. Donders Centre';
  content = 'timestamps from DMA log file';
elseif filetype_check_extension(filename, '.ttl') && filetype_check_header(filename, 'ttl')
  ftype = 'neuralynx_ttl';
  manufacturer = 'F.C. Donders Centre';
  content = 'Parallel_in from DMA log file';
elseif filetype_check_extension(filename, '.sdma') && isdir(filename)
  ftype = 'neuralynx_sdma';
  manufacturer = 'F.C. Donders Centre';
  content = 'split DMA log file';

  % known Plexon file types
elseif filetype_check_extension(filename, '.nex')  && filetype_check_header(filename, 'NEX1')
  ftype = 'plexon_nex';
  manufacturer = 'Plexon';
  content = 'electrophysiological data';
elseif filetype_check_extension(filename, '.plx')  && filetype_check_header(filename, 'PLEX')
  ftype = 'plexon_plx';
  manufacturer = 'Plexon';
  content = 'electrophysiological data';
elseif filetype_check_extension(filename, '.ddt')
  ftype = 'plexon_ddt';
  manufacturer = 'Plexon';
elseif isdir(filename) && most(filetype_check_extension({ls.name}, '.nex')) && most(filetype_check_header({ls.name}, 'NEX1'))
  % a directory containing multiple plexon NEX files
  ftype = 'plexon_ds';
  manufacturer = 'Plexon';
  content = 'electrophysiological data';

  % known Cambridge Electronic Design file types
elseif filetype_check_extension(filename, '.smr')
  ftype = 'ced_son';
  manufacturer = 'Cambridge Electronic Design';
  content = 'Spike2 SON filing system';

  % known BESA file types
elseif filetype_check_extension(filename, '.avr') && strcmp(ftype, 'unknown')
  ftype = 'besa_avr';  % FIXME, can also be EEProbe average EEG
  manufacturer = 'BESA';
  content = 'average EEG';
elseif filetype_check_extension(filename, '.elp')
  ftype = 'besa_elp';
  manufacturer = 'BESA';
  content = 'electrode positions';
elseif filetype_check_extension(filename, '.eps')
  ftype = 'besa_eps';
  manufacturer = 'BESA';
  content = 'digitizer information';
elseif filetype_check_extension(filename, '.sfp')
  ftype = 'besa_sfp';
  manufacturer = 'BESA';
  content = 'sensor positions';
elseif filetype_check_extension(filename, '.ela')
  ftype = 'besa_ela';
  manufacturer = 'BESA';
  content = 'sensor information';
elseif filetype_check_extension(filename, '.pdg')
  ftype = 'besa_pdg';
  manufacturer = 'BESA';
  content = 'paradigm file';
elseif filetype_check_extension(filename, '.tfc')
  ftype = 'besa_tfc';
  manufacturer = 'BESA';
  content = 'time frequency coherence';
elseif filetype_check_extension(filename, '.mul')
  ftype = 'besa_mul';
  manufacturer = 'BESA';
  content = 'multiplexed ascii format';
elseif filetype_check_extension(filename, '.gen')
  ftype = 'besa_gen';
  manufacturer = 'BESA';
  content = 'generic ascii format';
elseif filetype_check_header(filename, 'BESA_SA_IMAGE')
  ftype = 'besa_src';
  manufacturer = 'BESA';
  content = 'beamformer source reconstruction';
elseif filetype_check_extension(filename, '.swf') && filetype_check_header(filename, 'Npts=')
  ftype = 'besa_swf';
  manufacturer = 'BESA';
  content = 'beamformer source waveform';

  % files from Pascal Fries' PhD research at the MPI
elseif filetype_check_extension(filename, '.dap') && filetype_check_header(filename, char(1))
  ftype = 'mpi_dap';
  manufacturer = 'MPI Frankfurt';
  content = 'electrophysiological data';
elseif isdir(filename) && ~isempty(cell2mat(regexp({ls.name}, '.dap$')))
  ftype = 'mpi_ds';
  manufacturer = 'MPI Frankfurt';
  content = 'electrophysiological data';

  % known Curry V4 file types
elseif filetype_check_extension(filename, '.dap')
  ftype = 'curry_dap';   % FIXME, can also be MPI Frankfurt electrophysiological data
  manufacturer = 'Curry';
  content = 'data parameter file';
elseif filetype_check_extension(filename, '.dat')
  ftype = 'curry_dat';
  manufacturer = 'Curry';
  content = 'raw data file';
elseif filetype_check_extension(filename, '.rs4')
  ftype = 'curry_rs4';
  manufacturer = 'Curry';
  content = 'sensor geometry file';
elseif filetype_check_extension(filename, '.par')
  ftype = 'curry_par';
  manufacturer = 'Curry';
  content = 'data or image parameter file';
elseif filetype_check_extension(filename, '.bd0') || filetype_check_extension(filename, '.bd1') || filetype_check_extension(filename, '.bd2') || filetype_check_extension(filename, '.bd3') || filetype_check_extension(filename, '.bd4') || filetype_check_extension(filename, '.bd5') || filetype_check_extension(filename, '.bd6') || filetype_check_extension(filename, '.bd7') || filetype_check_extension(filename, '.bd8') || filetype_check_extension(filename, '.bd9')
  ftype = 'curry_bd';
  manufacturer = 'Curry';
  content = 'BEM description file';
elseif filetype_check_extension(filename, '.bt0') || filetype_check_extension(filename, '.bt1') || filetype_check_extension(filename, '.bt2') || filetype_check_extension(filename, '.bt3') || filetype_check_extension(filename, '.bt4') || filetype_check_extension(filename, '.bt5') || filetype_check_extension(filename, '.bt6') || filetype_check_extension(filename, '.bt7') || filetype_check_extension(filename, '.bt8') || filetype_check_extension(filename, '.bt9')
  ftype = 'curry_bt';
  manufacturer = 'Curry';
  content = 'BEM transfer matrix file';
elseif filetype_check_extension(filename, '.bm0') || filetype_check_extension(filename, '.bm1') || filetype_check_extension(filename, '.bm2') || filetype_check_extension(filename, '.bm3') || filetype_check_extension(filename, '.bm4') || filetype_check_extension(filename, '.bm5') || filetype_check_extension(filename, '.bm6') || filetype_check_extension(filename, '.bm7') || filetype_check_extension(filename, '.bm8') || filetype_check_extension(filename, '.bm9')
  ftype = 'curry_bm';
  manufacturer = 'Curry';
  content = 'BEM full matrix file';
elseif filetype_check_extension(filename, '.dig')
  ftype = 'curry_dig';
  manufacturer = 'Curry';
  content = 'digitizer file';

  % known Curry V2 file types
elseif filetype_check_extension(filename, '.sp0') || filetype_check_extension(filename, '.sp1') || filetype_check_extension(filename, '.sp2') || filetype_check_extension(filename, '.sp3') || filetype_check_extension(filename, '.sp4') || filetype_check_extension(filename, '.sp5') || filetype_check_extension(filename, '.sp6') || filetype_check_extension(filename, '.sp7') || filetype_check_extension(filename, '.sp8') || filetype_check_extension(filename, '.sp9')
  ftype = 'curry_sp';
  manufacturer = 'Curry';
  content = 'point list';
elseif filetype_check_extension(filename, '.s10') || filetype_check_extension(filename, '.s11') || filetype_check_extension(filename, '.s12') || filetype_check_extension(filename, '.s13') || filetype_check_extension(filename, '.s14') || filetype_check_extension(filename, '.s15') || filetype_check_extension(filename, '.s16') || filetype_check_extension(filename, '.s17') || filetype_check_extension(filename, '.s18') || filetype_check_extension(filename, '.s19') || filetype_check_extension(filename, '.s20') || filetype_check_extension(filename, '.s21') || filetype_check_extension(filename, '.s22') || filetype_check_extension(filename, '.s23') || filetype_check_extension(filename, '.s24') || filetype_check_extension(filename, '.s25') || filetype_check_extension(filename, '.s26') || filetype_check_extension(filename, '.s27') || filetype_check_extension(filename, '.s28') || filetype_check_extension(filename, '.s29') || filetype_check_extension(filename, '.s30') || filetype_check_extension(filename, '.s31') || filetype_check_extension(filename, '.s32') || filetype_check_extension(filename, '.s33') || filetype_check_extension(filename, '.s34') || filetype_check_extension(filename, '.s35') || filetype_check_extension(filename, '.s36') || filetype_check_extension(filename, '.s37') || filetype_check_extension(filename, '.s38') || filetype_check_extension(filename, '.s39')
  ftype = 'curry_s';
  manufacturer = 'Curry';
  content = 'triangle or tetraedra list';
elseif filetype_check_extension(filename, '.pom')
  ftype = 'curry_pom';
  manufacturer = 'Curry';
  content = 'anatomical localization file';
elseif filetype_check_extension(filename, '.res')
  ftype = 'curry_res';
  manufacturer = 'Curry';
  content = 'functional localization file';

  % known MBFYS file types
elseif filetype_check_extension(filename, '.tri')
  ftype = 'mbfys_tri';
  manufacturer = 'MBFYS';
  content = 'triangulated surface';
elseif filetype_check_extension(filename, '.ama') && filetype_check_header(filename, [0 0 0 10])
  ftype = 'mbfys_ama';
  manufacturer = 'MBFYS';
  content = 'BEM volume conduction model';

  % some other known file types
elseif length(filename>4) && exist([filename(1:(end-4)) '.mat'], 'file') && exist([filename(1:(end-4)) '.bin'], 'file')
  % this is a self-defined FCDC data format, consisting of two files
  % there is a matlab V6 file with the header and a binary file with the data (multiplexed, ieee-le, double)
  ftype = 'fcdc_matbin';
  manufacturer = 'F.C. Donders Centre';
  content = 'multiplexed ielectrophysiology data';
elseif filetype_check_extension(filename, '.lay')
  ftype = 'layout';
  manufacturer = 'Ole Jensen';
  content = 'layout of channels for plotting';
elseif filetype_check_extension(filename, '.dcm') || filetype_check_extension(filename, '.ima')
  ftype = 'dicom';
  manufacturer = 'Dicom';
  content = 'image data';
elseif filetype_check_extension(filename, '.trl')
  ftype = 'fcdc_trl';
  manufacturer = 'F.C.Donders';
  content = 'trial definitions';
elseif filetype_check_header(filename, [255 'BIOSEMI']) % filetype_check_extension(filename, '.bdf')
  ftype = 'biosemi_bdf';
  %   ftype = 'bham_bdf';
  manufacturer = 'Biosemi Data Format';
  content = 'electrophysiological data';
elseif filetype_check_extension(filename, '.edf')
  ftype = 'edf';
  manufacturer = 'European Data Format';
  content = 'electrophysiological data';
elseif filetype_check_extension(filename, '.mat') && filetype_check_header(filename, 'MATLAB')
  ftype = 'matlab';
  manufacturer = 'Matlab';
  content = 'Matlab binary data';
end

if strcmp(ftype, 'unknown')
  warning(sprintf('could not determine filetype of %s', filename));
end

% remember the details
detail.filetype     = ftype;
detail.manufacturer = manufacturer;
detail.content      = content;

if ~isempty(desired)
  % return a boolean value instead of a descriptive string
  ftype = strcmp(ftype, desired);
end

% SUBFUNCTION that helps in deciding whether a directory with files should
% be treated as a "dataset". This function returns a logical 1 (TRUE) if more 
% than half of the element of a vector are nonzero number or are logical 1 (TRUE).
function y = most(x);
x = x(find(~isnan(x(:))));
y = sum(x==0)<ceil(length(x)/2);
