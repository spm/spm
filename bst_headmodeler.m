function [varargout] = bst_headmodeler(varargin);
%BST_HEADMODELER - Solution to the MEG/EEG forward problem
% function [varargout] = bst_headmodeler(varargin);
% Authorized syntax:
% [G, OPTIONS] = bst_headmodeler(StudyFile, OPTIONS);
% [G, OPTIONS] = bst_headmodeler(OPTIONS);
% [OPTIONS] = bst_headmodeler;
%
% --------------------------------- INPUTS -------------------------------------
% INPUTS
%
%      [OPTIONS] = bst_headmodeler; Returns the default values for the OPTIONS
%      parameter structure
% 
%      StudyFile: the name of a BrainStorm study file. We suppose the
%      corresponding BrainStorm channel file  is available in the same folder
%      with the conventional file name (e.g., for a study file calles
%      meg_brainstormstudy.mat, BST_HEADMODELER expects to find the
%      corresponding channel file under the name meg_channel.mat). If the
%      channel file were not sitting in the same folder as  the study file,
%      OPTIONS.ChannelFile enforces computation of the forward model with the
%      channel information contained in OPTIONS.ChannelFile.
%      
%      If no further input arguments are specified, forward modeling is
%      completed with the default parameters specified below
%         
%      OPTIONS a structure where optional parameters may be specified using the
%      following fields  Note: if no options are specified, BST_HEADMODELER will
%      proceed to the computation of the foward problem on a 3D grid of source
%      locations that cover the entire head volume See
%      OPTIONS.VolumeSourceGridSpacing for default settings
%  
%      Important notice: there is no need to define all the following fields
%      when using the OPTIONS argument. The undefined field(s) will be assigned
%      default values.
%
%      *
%      * Fields Related to forward approach
%      *
%
%      .Method: is either a character string or a cell array of two strings that
%               specifies the kind of approach to be applied to the compuation
%               of the foward model In case Method is a cell array, it should
%               contain 2 strings, one to specifiy the method to be used for MEG
%               and the other for . If only a single string is specified, the
%               foward computation will be completed on the set of corresponding
%               CHANndx only (i.e. MEG or MEG) Available forward modeling
%               methods and corresponding authorized strings for Method:
%               - MEG
%                    'meg_sphere' (DEFAULT) : Spherical head model designed
%                       following the Sarvas analytical formulation (i.e.
%                       considering the true orientation 
%                       of the magnetic field sensors) (see OPTIONS.HeadCenter)
%                    'meg_os' : MEG overlapping sphere forward model
%                    'meg_bem' : Apply BEM computation (see OPTIONS.BEM for details)
%               - 
%                    'eeg_sphere' : Single-sphere forward modeling (see
%                       OPTIONS.HeadCenter, OPTIONS.Radii, OPTIONS.Conductivity)
%                    'eeg_3sphere' : EEG forward modeling with a set of 3
%                       concentric spheres (Scalp, Skull, Brain/CSF) (see
%                       OPTIONS.HeadCenter, OPTIONS.Radii, OPTIONS.Conductivity)
%                    'eeg_3sphereBerg' (DEFAULT) : Same as eeg_3sphere with
%                       correction for possible dipoles outside the sphere
%                    'eeg_os' : EEG overlapping sphere head model (see
%                       OPTIONS.HeadCenter, OPTIONS.Radii, OPTIONS.Conductivity)
%                    'eeg_bem' : Apply BEM computation (see OPIONS.BEM for details)
% 
%                Default is {'meg_sphere','eeg_3sphereBerg'};   
%
%      .HeadModelName : a character string that specifies the name of the
%                       headmodel represented by this file, e.g "Spherical", 
%                       "Overlapping Spheres", "Constant Collocation BEM", etc.
%                       Default is "Default", meaning it will include the the
%                       name(s) of the method(s) used in the MEG and/or EEG
%                       forward models
%
%      *
%      * Fields Related to function's I/O
%      *
%      
%      .HeadModelFile : Specifies the name of the head model file where to store
%                       the forward model.  If set to 'default', the default
%                       nomenclature for BrainStorm's head model file name is
%                       used and BST_HEADMODELER creates a file in StudyFile's
%                       folder.
%                       Default is empty. 
%      .ImageGridFile : Specifies the name of the file where to store the full
%                       cortical gain matrix file  If set to 'default', the
%                       default nomenclature for BrainStorm's head model file
%                       name is used and BST_HEADMODELER creates a file in
%                       StudyFile's folder.
%                       Default is empty.
%      .ImageGridBlockSize : Number of sources for which to compute the forward
%                            model at a time in a block computation routines
%                            (saves memory space). This option is relevant only
%                            when some forward modeling on cortical surface is
%                            requested (i.e. when .Cortex is specified)
%                            Default is 2000
%      .FileNamePrefix : A string that specifies the prefix for all file names
%                       (Channel, HeadModel, Gain Matrices) when .HeadModelFile
%                       is set to 'default' and .ChannelFile is empty.
%                        Default is 'bst_'
%      .Verbose     : Toggles verbose mode on/off;
%                     Default is 1 (on)
%    
%      *
%      * Fields Related to Head Geometry *
%      *
% 
%      .Scalp     : A structure specifying the Scalp surface envelope to serve
%                   for parameter adjustment of best-fitting sphere, with
%                   following fields:
%                   .FileName : A string specifying the name of the BrainStorm
%                               tessellation file containing the Scalp
%                               tessellation (default is 1);
%                   .iGrid    : An integer for the index of the Scalp surface in
%                               the Faces, Vertices and Comments cell arrays in
%                               the tessellation file
%                               Default is empty (Best-fitting sphere is
%                               computed from the sensor array).
%      .HeadCenter: a 3-element vector specifying the coordinates, in the
%                   sensors coordinate system, of the center of the spheres that
%                   might be used in the head model. 
%                   Default is estimated from the center of the best-fitting
%                   sphere to the sensor locations
%      .Radii     : a 3-element vector containing the radii of the single or 3
%                   concentric spheres, when needed;
%                   Order must be the following : [Rcsf, Routerskull, Rscalp];
%                   Default is estimated from the best-fitting sphere to the
%                   sensor locations and  OPTIONS.Radii is set to: Rscalp [.88
%                   .93 1]. Rscalp is estimated from the radius of the
%                   best-fitting sphere;
%      .Conductivity : a 3-element vector containing the values for the
%                      conductivity of the tissues in the following order:
%                      [Ccsf, Cskull, Cscalp];
%                      Default is set to [.33 .0042 .33];
%      .EEGRef    : the NAME (not index of the channel file) of the electrode
%                   that acts as the reference channel for the EEG. If data is
%                   referenced to instantaneous average (i.e. so called
%                   average-reference recording) value is 'AVERAGE REF'; 
%                   IMPORTANT NOTICE: When user calls bst_headmodeler with the
%                   .ChannelLoc option and .ChannelType = 'EEG'and wants the EEG
%                   reference to be e.g. channel 26, then .EEGRef should be set
%                   to 'EEG 26'                
%                   Default is 'AVERAGE REF'.
%      
%      .OS_ComputeParam : if 1, force computation of all sphere parameters when
%                         choosing a method based on either the MEG or EEG
%                         overlapping-sphere technique,  if 0 and when
%                         .HeadModelFile is specified, sphere parameters are
%                         loaded from the pre-computed HeadModel file.
%                         Default is 1.
%
%      .BEM       : Structure that specifies the necessary BEM parameters    
%                   .Interpolative : Flag indicating whether exact or
%                                    interpolative approach is used to compute
%                                    the forward solution using BEM. 
%                                    if set to 1, exact computation is run on a
%                                    set of points distributed wihtin the inner
%                                    head volume and any subsequent request for
%                                    a forward gain vector (e.g. during a
%                                    volumic source scan using RAP-MUSIC) is
%                                    computed using an interpolation of the
%                                    forward gain vectors of the closest volumic
%                                    grid points. This allows faster computation
%                                    of the BEM solution during source search.
%                                    if set to 0, exact computation is required
%                                    at every source location.
%                                    We recommend to set it to 0 (default) when
%                                    sources have fixed location, e.g.
%                                    constrained on the cortical surface.
%                   .EnvelopeNames : a cell array of strutures that specifies
%                                    the ORDERED tessellated surfaces to be
%                                    included in the BEM computation. 
%                                    .EnvelopeNames{k}.TessFile : A string for
%                                    the name of the tessellation file
%                                    containing the kth surface
%                                    .EnvelopeNames{k}.TessName : A string for
%                                    the name of the surface within the
%                                    tessellation file  This string should match
%                                    one of the Comment strings in the
%                                    tessellation file. The chosen surfaces must
%                                    be ordered starting by the the innermost
%                                    surface  (e.g. brain or inner skull
%                                    surface) and finishing with the outermost
%                                    layer (e.g. the scalp)
%                   .Basis         : set to either 'constant' or 'linear' (default)
%                   .Test          : set to either 'Galerkin' or 'Collocation' (default)
%                   .ISA           : Isolated-skull approach set to 0 or 1 (default is 1)
%                   .NVertMax      : Maximum number of vertices per envelope,
%                                    therefore leading to decimation of orginal
%                                    surfaces if necessary
%                                    (default is 1000)
%                   .ForceXferComputation: if set to 1, force recomputation of
%                                          existing transfer matrices in current
%                                          study folder (replace existing
%                                          files);
%                                          Default is 1
%
%      *
%      * Fields Related to Sensor Information *
%      *
%      
%      .ChannelFile : Specifies the name of the file containing the channel
%                     information (needs to be a BrainStorm channel file). If
%                     file does not exists and sensor information is provided in
%                     .ChannelLoc, a BrainStorm Channl file with name
%                     .ChannelFile is created
%                     Default is left blank as this information is extracted
%                     from the channel file associated to the chosen BrainStorm
%                     studyfile.
%      .Channel     : A full BrainStorm channel structure if no channel file is
%                     specified;
%                     Default is empty 
%      .ChannelType : A string specifying the type of channel in ChannelLoc. Can
%                     be either 'MEG' or 'EEG'. Note that the same channel type
%                     is assumed for every channel.
%                     Default is empty
%      .ChannelLoc  : Specifies the location of the channels where to compute
%                     the forward model Can be either a 3xNsens (for EEG or
%                     MEG-magnetometer) or 6xNsens matrix (for the
%                     MEG-gradiometer case). (for magnetometer or gradiometer
%                     MEG - channel weights are set to -1 and 1 for each
%                     magnetometer in the gradiometer respectively)  Note that
%                     in the MEG-gradiometer case, the 3 first (res. last) rows
%                     of .ChannelLoc stand for each of the magnetometers of the
%                     gradiometer set. In the case of a mixture of MEG magneto
%                     and MEG gradio-meters, .ChannelLoc needs to be a 6xNsens
%                     matrix where the last 3 rows are filled with NaN for
%                     MEG-magnetometers. If ones wants to handle both EEG and MEG
%                     sensors, please create a full ChannelFile and use the
%                     .ChannelFile option.
%                     Default is empty (Information extracted from the ChannelFile).
%      .ChannelOrient : Specifies the orientation of the channels where to
%                       compute the EEG or MEG forward model   Can be either a
%                       3xNsens (for EEG and magnetometer MEG) or 6xNsens (for
%                       gradiometer MEG) matrix or a cell array of such matrices
%                       (one cell per type of method selected)
%                       Default is empty (Information extracted from the
%                       ChannelFile or assume radial orientation when
%                       .ChannelLoc is filled).
%
%      *
%      * Fields Related to Source Models *
%      *
%      .SourceModel : A vector indicating the type of source models to be
%                     computed; The following code is enfoced:
%                     -1 : Compute the forward fields of Current Dipole sources
%                          (available for all forward approaches)
%                      1 : 1st-order Current Multipole Sources 
%                          (available for sphere-based MEG approaches only)
%                     User can set OPTIONS.SourceModel to e.g., [-1 1] to
%                     compute forward models from both source models.
%                     Default is -1 
%
%
%      *
%      * Fields Related to Source Localization *
%      *
%      .Cortex                  :  A structure specifying the Cortex surface
%                                  envelope to serve as an image support with
%                                  following fields.
%                                  .FileName : A string specifying the name of
%                                              the BrainStorm tessellation file
%                                              containing the Cortex
%                                              tessellation;
%                                  .iGrid    : An integer for the index of the
%                                              Cortex surface in the Faces,
%                                              Vertices and Comments cell arrays
%                                              in the tessellation file (default
%                                              is 1)
%                                              Default is empty.
%      .GridLoc                 :  A 3xNsources matrix that contains the
%                                  locations of the sources at which the forward
%                                  model will be computed
%                                  Default is empty (Information taken from
%                                  OPTIONS.Cortex or OPTIONS.VolumeSourceGrid);
%      .GridOrient              :  a 3xNsources matrix that forces the source
%                                  orientation at every vertex of the .ImageGrid
%                                  cortical surface;
%                                  Defaults is empty; this information being
%                                  extracted from the corresponding tessellated
%                                  surface. 
%      .ApplyGridOrient         : if set to 1, force computation of the forward
%                                 fields by considering the local orientation of
%                                 the cortical surface;
%                                 If set to 0, a set of 3 orthogonal dipoles are
%                                 considered at each vertex location on the
%                                 tessellated surface.      
%                                 Default is 1.
%      .VolumeSourceGrid        : if set to 1, a 3D source grid is designed to
%                                 fit inside the head volume and will serve as a
%                                 source space for scannig techniques such as
%                                 RAP-MUSIC;
%                                 if set to 0, this grid will be computed at the
%                                 first call of e.g. RAP-MUSIC);
%                                 Default is 1 
%      .VolumeSourceGridSpacing : Spacing in centimeters between two consecutive
%                                 sources in the 3D source grid described above;
%                                 Default is 2 cm.
%      .VolumeSourceGridLoc     : a 3xN matrix specifying the locations of the
%                                 grid points that will be used to design the
%                                 volumic search grid (see .VolumicSourceGrid)
%                                 Default is empty (locations are estimated
%                                 automatically to cover the estimated inner
%                                 head volume)
%      .SourceLoc               : a 3xNsources matrix that contains the
%                                 locations of the sources at which the forward
%                                 model will be computed
%                                 Default is empty (Information taken from
%                                 OPTIONS.ImageGrid or
%                                 OPTIONS.VolumeSourceGrid);
%      .SourceOrient            : a 3xNsources matrix that contains the
%                                 orientations of the sources at which the
%                                 forward model will be computed
%                                 Default is empty.
%
%
% --------------------------------- OUTPUTS ------------------------------------
% OUTPUT
%      G        if the number of sources (Nsources) if less than
%               .ImageGridBlockSize then G is a gain matrix of dimension
%               Nsensors x Nsources: Each column of G is the forward field
%               created by a dipolar source of unit amplitude. Otherwise, G is
%               the name of the binary file containing the gain matrix. This
%               file can be read using the READ_GAIN function.
%      OPTIONS  Returns the OPTIONS structure with updated fields following the
%               call to BST_HEADMODELER.  Can be useful to obtain a full
%               BrainStorm Channel structure when only the .ChannelLoc and
%               possibly .ChannelOrient fields were provided.

%<autobegin> ---------------------- 27-Jun-2005 10:43:31 -----------------------
% ------ Automatically Generated Comments Block Using AUTO_COMMENTS_PRE7 -------
%
% CATEGORY: Forward Modeling
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\bem_gain.m
%   toolbox\bem_xfer.m
%   toolbox\berg.m
%   toolbox\bst_message_window.m
%   toolbox\colnorm.m
%   toolbox\get_channel.m
%   toolbox\get_user_directory.m
%   toolbox\good_channel.m
%   toolbox\gridmaker.m
%   toolbox\inorcol.m
%   toolbox\norlig.m
%   toolbox\overlapping_sphere.m
%   toolbox\rownorm.m
%   toolbox\save_fieldnames.m
%   toolbox\source_grids.m
%   toolbox\view_surface.m
%
% Subfunctions in this file, in order of occurrence in file:
%   BEMGaingridFname = bem_GainGrid(DataType,OPTIONS,BEMChanNdx)
%   g = gterm_constant(r,rq)
%
% At Check-in: $Author: Silvin $  $Revision: 68 $  $Date: 12/15/05 4:14a $
%
% This software is part of BrainStorm Toolbox Version 27-June-2005  
% 
% Principal Investigators and Developers:
% ** Richard M. Leahy, PhD, Signal & Image Processing Institute,
%    University of Southern California, Los Angeles, CA
% ** John C. Mosher, PhD, Biophysics Group,
%    Los Alamos National Laboratory, Los Alamos, NM
% ** Sylvain Baillet, PhD, Cognitive Neuroscience & Brain Imaging Laboratory,
%    CNRS, Hopital de la Salpetriere, Paris, France
% 
% See BrainStorm website at http://neuroimage.usc.edu for further information.
% 
% Copyright (c) 2005 BrainStorm by the University of Southern California
% This software distributed  under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html .
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%<autoend> ------------------------ 27-Jun-2005 10:43:31 -----------------------



% /---Script Author--------------------------------------\
% |                                                      |
% | *** Sylvain Baillet Ph.D.                            |
% | Cognitive Neuroscience & Brain Imaging Laboratory    |
% | CNRS UPR640 - LENA                                   | 
% | Hopital de la Salpetriere, Paris, France             |
% | sylvain.baillet@chups.jussieu.fr                     |
% |                                                      |
% \------------------------------------------------------/
%  
% Date of creation: March 2002


% ----------------------------- Script History ---------------------------------
%
% SB 07-Aug-2002  Fixed GridLoc saving issues - now GridLoc
%                 is a cell array of length one which indicates
%                 the name of the tessellation file that was
%                 used as an imaging support.
% SB 03-Sep-2002  Fixed minor bugs when function is called from command line
%                 Now bst_headmodeler can be called even when Brainstorm is not running
% SB 05-Sep-2002  BEM can now be called from command line as well. 
%                 Updated script header
% SB 06-Sep-2002  Updated EEG BEM : added option for arbitrary reference of the
%                 potentials (takes now OPTIONS.EEGRef into account).
% SB 18-Nov-2002  Updated EEG reference management.
% SB 21-Nov-2002  Minor fixes for command line calls.
% SB 21-Nov-2002bis Fixed bug from the get(*,'VertexNormals') command that
%                 sometimes returns [0 0 0] at some vertex                 
% SB 23-Dec-2002  Correct HeadModelFile name is passed to OPTIONS output argument
% JCM 27-May-2004 Cleaning comments
%                 Found Matlab bug and temp solution.
%                 Matlab 6.5.0 bug, Technical Solution Number: 1-19NOK
%                 http://www.mathworks.com/support/solutions/data/1-19NOK.html?solution=1-19NOK
%                 Must REWIND first before fseek.
%                 Apparently fixed in 6.5.1
%                 Example:
%                 status = fseek(fdest,0,'bof');
%                 status = fseek(fdest,offset,'bof');
% SB  08-Mar-2005 Fixed bug in BEM computation 
% SB  02-Feb-2005 Fixed bug that occured when ChannelLoc field is used to
%                 passed EEG channel locations (thanks Jeremie Mattout at the FIL)
% ----------------------------- Script History ---------------------------------

% Default options settings--------------------------------------------------------------------------------------------------
DefaultMethod = {'meg_sphere','eeg_3sphereBerg'};

ReducePatchScalpNVerts = 500; % Number of vertices the scalp tessellation will be reduced to in OS computations
BEM_defaults = struct(...
    'Basis','linear',...
    'Test','galerkin',...
    'Interpolative',0,...
    'ISA',1,...
    'NVertMax',1000,...
    'ForceXferComputation', 1, ...
    'checksurf',0);

Def_OPTIONS = struct(...
    'ApplyGridOrient',1,...
    'BEM', BEM_defaults,...
    'Channel', [],...
    'ChannelFile', '',...
    'ChannelLoc', '',...
    'ChannelOrient', '',...
    'ChannelType', '',...
    'Conductivity', [.33 .0042 .33],...
    'Cortex',[],...
    'EEGRef','',...
    'HeadCenter',[],...
    'HeadModelFile', '',...
    'HeadModelName','Default',...
    'ImageGridBlockSize', 2000,...
    'ImageGridFile', '',...
    'GridOrient',[],...
    'Method', {DefaultMethod},...
    'OS_ComputeParam', 1,...    
    'PrefixFileName','bst_',...
    'Radii', [],...
    'Scalp',[],...    
    'SourceLoc',[],...
    'SourceModel', [-1],...
    'SourceOrient',[],...
    'StudyFile','',...
    'TessellationFile','',...
    'VolumeSourceGrid',1,...
    'VolumeSourceGridSpacing', 2,...
    'VolumeSourceGridLoc', [],...
    'Verbose', 1 ...
    );

SourceOrderString = {'Current Dipole',...
        'Unvalid Order (was Magnetic Dipole Moment)',...
        'Current Multipole Expansion',...
        'Current Dipole Pairs',...
        'Unvalid Order (was Magnetic Dipole Moment PAIRS)',...    
        'Current Multipole Expansion Pairs'};

if nargin == 0
    if nargout > 1
        varargout{1} = Def_OPTIONS;
        varargout{2} = Def_OPTIONS;
    else
        varargout{1} = Def_OPTIONS;
    end
    return
elseif nargin == 1
    if ischar(varargin{1}) % User enters only the studyfile name
        StudyFile = varargin{1};
        OPTIONS = Def_OPTIONS;
    elseif isstruct(varargin{1}) % User enters only the OPTION field
        OPTIONS = varargin{1};
    else
        errordlg('Uncorrect input parameter type, please check the help file'), varargout = cell(nargout,1);
        return, 
    end
elseif nargin == 2 % No options were specified, apply default
    StudyFile = varargin{1};
    OPTIONS = varargin{2};
else
    errordlg('Wrong number of arguments when calling head modeler')
    varargout = cell(nargout,1);
    return
end

% Check field names of passed OPTIONS and fill missing ones with default values
DefFieldNames = fieldnames(Def_OPTIONS);
for k = 1:length(DefFieldNames)
    if ~isfield(OPTIONS,DefFieldNames{k}) | strcmp(DefFieldNames{k},'BEM')
        if ~isfield(OPTIONS,DefFieldNames{k})
            
            OPTIONS = setfield(OPTIONS,DefFieldNames{k},getfield(Def_OPTIONS,DefFieldNames{k}));
            
        elseif strcmp(DefFieldNames{k},'BEM')
            
            BEM_DefFieldNames = fieldnames(BEM_defaults);
            for kk = 1:length(BEM_DefFieldNames)
                if ~isfield(OPTIONS.BEM,BEM_DefFieldNames{kk})
                    OPTIONS.BEM = setfield(OPTIONS.BEM,BEM_DefFieldNames{kk},getfield(BEM_defaults,BEM_DefFieldNames{kk}));
                end
            end
        end
    end
end

if isempty(OPTIONS.Conductivity)
    OPTIONS.Conductivity = Def_OPTIONS.Conductivity;
end

clear Def_OPTIONS


if isempty(OPTIONS.HeadModelFile) & ~isempty(OPTIONS.ImageGridFile)
    % Force creation of a headmodel file
    OPTIONS.HeadModelFile = 'default';
end

OPTIONS.HeadModelFileOld = OPTIONS.HeadModelFile;

% What type of forward model (MEG and/or EEG) ?
DataType.MEG = strmatch('meg',OPTIONS.Method); % Rank of the respective forward method in cell array Method (ie could be Method = {'meg_bem','eeg_sphere'} or vice-versa)
DataType.EEG = strmatch('eeg',OPTIONS.Method);

if ~iscell(OPTIONS.Method)
    OPTIONS.Method = {OPTIONS.Method};
end

MegMethod = [];
if ~isempty(DataType.MEG)
    MegMethod = OPTIONS.Method{DataType.MEG}; % String indicating the forward method selected for MEG (res. EEG)
end
EegMethod = [];
if ~isempty(DataType.EEG)
    EegMethod = OPTIONS.Method{DataType.EEG};
end


% Check inputs integrity

%Source models
if ~isempty(find(OPTIONS.SourceModel == 0)) | ~isempty(find(abs(OPTIONS.SourceModel) > 1)) % unValid source models
    if ~isempty(DataType.MEG)
        errordlg('Valid source model orders for MEG are: -1 (Current Dipole) and 1 (fist-order Current Multipole Expansion)')
    end
    if ~isempty(DataType.EEG)
        errordlg('Valid source model order for EEG is: -1 (Current Dipole) ')
    end
    varargout = cell(nargout,1);
    return
end

% Source locations
if isempty(OPTIONS.SourceLoc) & ~OPTIONS.VolumeSourceGrid & isempty(OPTIONS.Cortex) % No source locations were specified
    errordlg('No source locations are specified. Please fill either one of the following fields of OPTIONS: .SourceLoc / .VolumeSourceGrid / .Cortex')
    varargout = cell(nargout,1);
    return
end


%--------------------------------------------------------------------------------------------------------------------------------------------
%
% HEAD MODELING BEGINS
%
%--------------------------------------------------------------------------------------------------------------------------------------------

if OPTIONS.Verbose, 
    
    clockk = fix(clock);
    time0 = clock;

    bst_message_window({...
            ' ',...
            '__________________________________________',...
            sprintf(...
            'Head modeling begins (%s - %dH %dmin %ds)',date,clockk(4:6))...
        })
    
end

User = get_user_directory;
if isempty(User) % Function is called from command line: BrainStorm is not running 
    % Use default folders
    User.SUBJECTS = pwd;
    User.STUDIES = pwd;
end

try
    cd(User.STUDIES)
catch
end



% Get all Subject and Study information------------------------------------------------------------------------------------------------------
if exist('StudyFile','var')
    if OPTIONS.Verbose, bst_message_window('Loading Study and Subject Information...'), end
    try 
        load(StudyFile) % Load study information
        [StudyPath,tmp,ext] = fileparts(StudyFile); clear tmp    
    catch 
        errordlg(['Could not find ',StudyFile,' in current Matlab path'])
        varargout = cell(nargout,1);
        return
    end
    
    if isempty(OPTIONS.StudyFile)
        OPTIONS.StudyFile = StudyFile;
    end
    
    SubjectFile = BrainStormSubject; clear BrainStormSubject
    
    OPTIONS.rooot = findstr(StudyFile,'brainstormstudy.mat');
    if isempty(OPTIONS.rooot)
        errordlg('Study file name should be of the form ''*brainstormstudy.mat''')
        varargout = cell(nargout,1);
        return
    end
    OPTIONS.rooot = strrep(StudyFile,'brainstormstudy.mat','');
    
    Study = load(StudyFile);
    %User = get_user_directory;
    try
        cd(User.SUBJECTS)
        Subject = load(Study.BrainStormSubject);
    catch 
        errordlg(sprintf('Please make sure the subject''s information in %s is available on this plateform',Study.BrainStormSubject)); return
    end
    
    if isempty(OPTIONS.TessellationFile) % No Tessellation file was passed to the headmodeler
        if isfield(Subject,'Tesselation')
            OPTIONS.TessellationFile = Subject.Tesselation; % Take default (OBSOLETE as for BsT MMII because we now consider all tessellation files in Subject's folder)
        else
            OPTIONS.TessellationFile = '';
            %[OPTIONS.TessellationFile,DataPopup, Leader] = find_brainstorm_files('tess',fullfile(Users.SUBJECTS,OPTIONS.subjectpath));
        end
    end
end

% Get Channel Information ------------------------------------------------------------------------------------------------------------------
if ~isempty(OPTIONS.Channel)
    Channel = OPTIONS.Channel;
else
    if isempty(OPTIONS.ChannelFile) & isempty(OPTIONS.ChannelLoc) % Load Channel file in current folder
        
        [Channel, ChannelFile] = get_channel(fullfile(User.STUDIES,OPTIONS.StudyFile));
        OPTIONS.ChannelFile = fullfile(fileparts(fullfile(User.STUDIES,OPTIONS.StudyFile)),ChannelFile);
       
       if isempty(ChannelFile)
           errordlg(sprintf('Channel file %s is missing. Please have it available it the current study folder.',OPTIONS.ChannelFile), 'Missing Channel File')
           return
       end
        
    elseif isempty(OPTIONS.ChannelLoc) & exist(OPTIONS.ChannelFile,'file') % If no specific channel locations are given and channel file exists, load the proper channel file
        
        OPTIONS.rooot = strrep(lower(OPTIONS.ChannelFile),'channel.mat','');
        try
            load(OPTIONS.ChannelFile)
        catch
            cd(User.STUDIES)
            load(OPTIONS.ChannelFile)
        end
        
    else % Create a dummy Channel structure  with Channel Locations (res. Orientations) specified in OPTIONS.ChannelLoc (res .ChannelOrient)
        
        if OPTIONS.Verbose, bst_message_window('Creating the channel structure from information in the OPTIONS fields. . .'), end
        
        % Get Channel Locations
        nchan = size(OPTIONS.ChannelLoc,2); % Number of channels
        Channel = struct('Loc',[],'Orient',[],'Comment','','Weight',[],'Type','','Name','');
        Channel(1:nchan) = deal(Channel);
        
        ChanType = upper(OPTIONS.ChannelType);
        if isempty(ChanType)
            errordlg('Please specify a channel type (i.e. MEG or EEG) in the ChannelType field of OPTIONS'), 
            bst_message_window('Please specify a channel type (i.e. MEG or EEG) in the ChannelType field of OPTIONS'), 
            varargout = cell(nargout,1);
            return
        end
        [Channel(:).Type] = deal(ChanType);
        if size(OPTIONS.ChannelLoc,1) == 6 % MEG Gradiometers or mixture gradio/magneto meters
            OPTIONS.ChannelLoc = reshape(OPTIONS.ChannelLoc,3,nchan*2);
            iGradFlag = 1; % Flag - gradiometer-type sensor set
        else
            iGradFlag = 0; % Flag - EEG or magnetometer-type sensor set
        end
        ichan = 1;
        for k = 1:nchan
            if iGradFlag
                Channel(k).Loc = OPTIONS.ChannelLoc(:,ichan:ichan+1);
                ichan = ichan+2;
            else
                if strcmp(ChanType,'MEG')
                    Channel(k).Loc = OPTIONS.ChannelLoc(:,ichan);
                %elseif strcmp(ChanType,'EEG') % Artificially add a dummy column full of zeros to each Channel(k).Loc 
                    %Channel(k).Loc = [OPTIONS.ChannelLoc(:,ichan) [0 0 0]'];;
                elseif strcmp(ChanType,'EEG') % Artificially add a dummy column full of zeros to each Channel(k).Loc 
                    Channel(k).Loc = [OPTIONS.ChannelLoc(:,ichan)];;
                end
                ichan = ichan+1;
            end
            Channel(k).Name = sprintf('%s %d',ChanType,k);
            Channel(k).Comment = int2str(k);
        end
        clear ichan k
        
        % Get Channel Orientations
        if isempty(OPTIONS.ChannelOrient) & strcmp(ChanType,'MEG') % No channel orientation were specified: use radial sensors
            if OPTIONS.Verbose, bst_message_window('Assign radial orientation to all Channels. . .'), end
            if isempty(OPTIONS.HeadCenter)
                if iGradFlag
                    Vertices = OPTIONS.ChannelLoc(:,1:2:end)';
                else
                    Vertices = OPTIONS.ChannelLoc';
                end
                
                nscalp = size(Vertices,1);
                if nscalp > 500 % 500 points is more than enough to compute scalp's best fitting sphere
                    Vertices = Vertices(unique(round(linspace(1,nscalp,500))),:); 
                    nscalp = size(Vertices,1); 
                end
                nmes = size(Vertices,1);
                
                % Run parameters fit --------------------------------------------------------------------------------------------------------------
                mass = mean(Vertices); % center of mass of the scalp vertex locations   
                R0 = mean(norlig(Vertices - ones(nscalp,1)*mass)); % Average distance between the center of mass and the scalp points
                vec0 = [mass,R0];
                
                [minn,brp] = fminsearch('dist_sph',vec0,[],Vertices);
                OPTIONS.HeadCenter = minn(1:end-1); 
                if isempty(OPTIONS.Radii)
                    OPTIONS.Radii = minn(end);
                    OPTIONS.Radii = minn(end)*[1/1.14 1/1.08 1];
                end
                
                if OPTIONS.Verbose, bst_message_window({...
                            sprintf('Center of the Sphere : %3.1f %3.1f %3.1f (cm)',100*OPTIONS.HeadCenter),...
                            sprintf('Radius : %3.1f (cm)',100*OPTIONS.Radii(3)),...
                            '-> DONE',' '})
                end
                
                clear minn brp mass R0 Vertices vec0
            end
            
            tmp = [Channel.Loc] - repmat(OPTIONS.HeadCenter',1,nchan*(1+(iGradFlag==1)));
            OPTIONS.ChannelOrient = tmp*inorcol(tmp); % Radial orientation for every channel
            clear tmp
        elseif ~isempty(OPTIONS.ChannelOrient) & strcmp(ChanType,'MEG') & iGradFlag
            OPTIONS.ChannelOrient = reshape(OPTIONS.ChannelOrient,3,nchan*2);
        end
        
        if strcmp(ChanType,'MEG')
            for k=1:nchan
                Channel(k).Orient = OPTIONS.ChannelOrient(:,((1+(iGradFlag==1))*k-1):(1+(iGradFlag==1))*k);
            end
        end
        
        clear k
        
        % Define Weights
        if iGradFlag
            [Channel(:).Weight] = deal([1 -1]);
        else
            [Channel(:).Weight] = deal([1]);
        end
        
        if strcmp(ChanType,'MEG') % Force no reference channels
            Channel(1).irefsens = [];
        end
        
        % New Channel stbemructure completed: save as a new channel file
        if ~isempty(OPTIONS.ChannelFile)
            %OPTIONS.ChannelFile = [OPTIONS.rooot,'channel.mat'];
            save(OPTIONS.ChannelFile,'Channel')
            if OPTIONS.Verbose, bst_message_window({...
                        sprintf('Channel file created in : %s',OPTIONS.ChannelFile),...
                        '-> DONE',' '}), end
        end
        
    end
    
    %load(OPTIONS.ChannelFile)
    OPTIONS.Channel = Channel;
end

% Find MEG and EEG channel indices
MEGndx = good_channel(Channel,[],'MEG');
EEGndx = good_channel(Channel,[],'EEG');
EEGREFndx = good_channel(Channel,[],'EEG REF');

MEG = ~isempty(MEGndx) & ~isempty(DataType.MEG); %  == 1 if MEG is requested and is available 
EEG = ~isempty(EEGndx) & ~isempty(DataType.EEG); 

if ~EEG & ~MEG
    errordlg('Please check that data (MEG or EEG) and channel types are compatible.')
    return
end

% Test whether CME models are requested for EEG
if ~isempty(find(OPTIONS.SourceModel == 1)) & EEG
    if MEG
        if OPTIONS.Verbose
            bst_message_window('wrap',...
                {'',...
                    'CME model not available for EEG - Skipping. . .',...
                    'Keeping it for MEG forward modeling only',...
                    ''...
                })
        end
    else % EEG only is requested
        errordlg(...
            {'CME model not available for EEG - Skipping. . .',...
                'Keeping it for MEG forward modeling only'...
            })
        return
    end
    
end

% Detect EEG reference - if none is specified in the .Comment or Type fields, 
% and if none was passed through OPTIONS.EEGRef 
% -> Apply average-referencing of the potentials by default.

if EEG 
    if isempty(OPTIONS.EEGRef)
        
        % EEG Reference Channel
        EEGREFndx = good_channel(Channel,[],'EEG REF');
        if isempty(EEGREFndx) % Average EEG reference anyway
            [Channel(EEGndx).Comment] = deal('AVERAGE REF');
        end
        
    else % EEG Reference is specified
        
        switch(OPTIONS.EEGRef)
        case 'AVERAGE REF'
            [Channel(:).Comment] = deal('AVERAGE REF');
           
        otherwise
           
            EEGREFndx = strmatch(OPTIONS.EEGRef,char(Channel(:).Name));
            if isempty(EEGREFndx)
                errordlg(sprintf(...
                    'No channel named ''%s'' was found amongst available EEG channels. Cannot use it as a reference for EEG.',OPTIONS.EEGRef...
                    ))
                return
            end
            
        end
        
    end
    
    %     if isempty(EEGREFndx) & ~MEG % No reference electrode was defined for the EEG stop if no MEG to compute...
    %         errordlg('Please make sure you have defined a reference for the EEG')
    %         return
    %     elseif isempty(EEGREFndx) & MEG %  Skip EEG forward modeling and proceed to MEG
    %         EEG = 0;
    %         if OPTIONS.Verbose
    %             bst_message_window({...
    %                     'No reference was defined for the EEG',...
    %                     'Skipping EEG modeling and completing MEG only...',...
    %                 })
    %         end
    %     end
    
end

if OPTIONS.Verbose, 
    if EEG & MEG
        if ~isempty(EEGREFndx)
            bst_message_window({'Channel count for current computation:',[int2str(length(MEGndx)), ' MEG Channels'],...
                    sprintf('%d EEG Channels with Reference: %s',length(EEGndx),Channel(EEGREFndx).Name),' '})
        else
            bst_message_window({'Channel count for current computation:',[int2str(length(MEGndx)), ' MEG Channels'],...
                    sprintf('%d EEG Channels with Average Reference',length(EEGndx)),' '})
        end
        
    elseif MEG
        bst_message_window({'Channel count for current computation:',[int2str(length(MEGndx)), ' MEG Channels'],...
                ' '})
    elseif EEG 
        if ~isempty(EEGREFndx)
            bst_message_window({'Channel count for current computation:',...
                    sprintf('%d EEG Channels with Reference: %s',length(EEGndx),Channel(EEGREFndx).Name),' '})
        else 
            bst_message_window({'Channel count for current computation:',...
                    sprintf('%d EEG Channels with Average Reference',length(EEGndx)),' '})
        end
        
    end
end

if MEG & isempty(MEGndx) % User wants to compute MEG but no MEG data is available
    errordlg('Sorry - No MEG data is available'), return
end

if EEG & isempty(EEGndx) % User wants to compute EEG but no EEG data is available
    errordlg('Sorry - No EEG data is available'), return
end


% Computation of parameters of the best-fitting sphere --------------------------------------------------------------------------------------------------------------
if length(findstr('bem',[OPTIONS.Method{:}])) ~= length(OPTIONS.Method)% Only if sphere-based head model is requested in any modality (MEG or EEG)
    
    % Best-fitting sphere parameters --------------------------------------------------------------------------------------------------------------
    if isempty(OPTIONS.HeadCenter) | isempty(OPTIONS.Radii)
        
        if OPTIONS.Verbose, bst_message_window('Estimating Center of the Head. . .'), end
        
        if isempty(OPTIONS.Scalp) % Best-fitting sphere is derived from the sensor array
            %             ans = questdlg('No scalp surface was selected - Do you want to use a spherical approximation of the head derived from the sensor locations instead ?',...
            %                 '','Yes','No','Yes');
            if EEG & length(EEGndx) > 9 % If EEG is available but a reasonable number of electrodes, use it to compute the sphere parameters
                ndx = [EEGndx]; % Compute Sphere parameters from EEG sensors only
            else
                ndx = [MEGndx];
            end
            
            Vertices = zeros(length(ndx),3);
            for k=1:length(ndx)
                Vertices(k,:) = Channel(ndx(k)).Loc(:,1)';
            end
            
        else % Best-fitting sphere is computed from the set of vertices of the scalp surface enveloppe
            
            try
                load(fullfile(User.SUBJECTS,OPTIONS.Scalp.FileName),'Vertices')
            catch
                load(OPTIONS.Scalp.FileName,'Vertices')
            end
            
            if ~isfield(OPTIONS.Scalp,'iGrid') % Apply Default
                OPTIONS.Scalp.iGrid = 1;
            end
            
            try
                Vertices = Vertices{OPTIONS.Scalp.iGrid}';
            catch
                errordlg(sprintf(...
                    'Tessellation file %s does not contain %d enveloppes',...
                    OPTIONS.Scalp.FileName,OPTIONS.Scalp.iGrid))
            end
        end
        
        nscalp = size(Vertices,1);
        nmes = size(Vertices,1);
        
        % Run parameters fit --------------------------------------------------------------------------------------------------------------
        mass = mean(Vertices); % center of mass of the scalp vertex locations   
        R0 = mean(norlig(Vertices - ones(nscalp,1)*mass)); % Average distance between the center of mass and the scalp points
        vec0 = [mass,R0];
        
        [SphereParams,brp] = fminsearch('dist_sph',vec0,[],Vertices);
        
        
        if OPTIONS.Verbose
            bst_message_window({...
                    sprintf('Center of the Sphere : %3.1f %3.1f %3.1f (cm)',100*SphereParams(1:end-1)'),...
                    sprintf('Radius : %3.1f (cm)',100*SphereParams(end)),...
                    '-> DONE',' '})
            
        end
        
    end % no head center and sphere radii were specified by user
    
    % Assign default values for sphere parameters 
    % if none were specified before
    if isempty(OPTIONS.Radii)
        OPTIONS.Radii = SphereParams(end)*[1/1.14 1/1.08 1];
    end
    
    if isempty(OPTIONS.HeadCenter)
        OPTIONS.HeadCenter = SphereParams(1:end-1)'; 
    end
    
    if isempty(OPTIONS.Conductivity)
        OPTIONS.Conductivity = [.33 .0042 .33];
    end
    
end % Use BEM Approach, so proceed to the selection of the envelopes

% ---------------------------------------------------------------------------------------------------------
%Create HeadModel Param structure
Param(1:length(Channel)) = deal(struct('Center',[],'Radii',[],'Conductivity',[],'Berg',[]));

if ~isempty(findstr('os',[OPTIONS.Method{:}])) % Overlapping-Sphere EEG or MEG is requested: load the scalp's tessellation
    
    if isempty(OPTIONS.Scalp)
        errordlg('Please specify a subject tessellation file for Scalp enveloppe in OPTIONS.Scalp when using Overlapping-Sphere forward approach');
        return
    end
    
    % load(fullfile(User.SUBJECTS,BrainStorm.SubjectTess),'Faces');
    try
        load(fullfile(User.SUBJECTS,OPTIONS.Scalp.FileName),'Faces','Vertices');
    catch
        load(OPTIONS.Scalp.FileName,'Faces','Vertices');
    end

    Faces = Faces{OPTIONS.Scalp.iGrid};
    Vertices = Vertices{OPTIONS.Scalp.iGrid}';
    
    if size(Vertices,1) > ReducePatchScalpNVerts % Reducepatch the scalp tessellation for fastest OS computation
        nfv = reducepatch(Faces,Vertices,2*ReducePatchScalpNVerts);
        if OPTIONS.Verbose, bst_message_window({...
                    sprintf('Decimated scalp tessellation from %d to %d vertices.',size(Vertices,1),size(nfv.vertices,1)),...
                    ' '});
        end
        clear Faces Vertices
        SubjectFV.vertices = nfv.vertices';
        SubjectFV.faces = nfv.faces; 
        clear nfv
    else
        SubjectFV.vertices = Vertices; % Available from the computation of the head center above
        clear Vertices
        SubjectFV.faces = Faces; clear Faces
    end
    
    if length(findstr('os',[OPTIONS.Method{:}])) == 2 % OS approach requested fro both MEG and EEG
        ndx = [MEGndx,EEGndx];    
        
    else
        if MEG
            if strcmpi('meg_os',OPTIONS.Method{DataType.MEG}) % OS for MEG only
                ndx = MEGndx;
            end
        end
        
        if EEG 
            if strcmpi('eeg_os',OPTIONS.Method{DataType.EEG}) % OS for EEG only
                ndx = [EEGndx, EEGREFndx];
            end
        end
    end
    if OPTIONS.Verbose
        bst_message_window({...
                ' ', 'Computing Overlapping-sphere model. . .'}) 
    end
    
    if isempty(OPTIONS.HeadModelFile) | OPTIONS.OS_ComputeParam
        
        % Compute all spheres parameters
        %------------------------------------------------------------------
        Sphere = overlapping_sphere(Channel(ndx),SubjectFV,OPTIONS.Verbose,OPTIONS.Verbose);
        if OPTIONS.Verbose
            bst_message_window({...
                    'Computing Overlapping-sphere model -> DONE',' '}) 
        end
        [Param(ndx).Center] = deal(Sphere.Center);
        [Param(ndx).Radii]  = deal(Sphere.Radius);
        [Param(ndx).Conductivity] = deal(OPTIONS.Conductivity);
        
    elseif exist(OPTIONS.HeadModelFile,'file')
        
        load(OPTIONS.HeadModelFile,'Param') % or use precomputed
        if OPTIONS.Verbose
            bst_message_window({...
                    sprintf('Sphere parameters loaded from : %s -> DONE',OPTIONS.HeadModelFile),' '}) 
        end
        
    else
    
        errordlg(sprintf('Headmodel file %s does not exist in Matlab''s search path',OPTIONS.HeadModelFile))
        return
        
    end
    
end

% Saving HeadModel's Param structure in the headmodel file
if strcmp(lower(OPTIONS.HeadModelFile),'default')
    if ~isfield(OPTIONS,'rooot')
        OPTIONS.rooot = OPTIONS.PrefixFileName;
    end
    OPTIONS.HeadModelFile = [OPTIONS.rooot,'headmodel.mat'];
    ifile = 0; 
    while exist(OPTIONS.HeadModelFile,'file') % Do not overwrite headmodel files
        ifile = ifile + 1;
        OPTIONS.HeadModelFile = [OPTIONS.rooot,'headmodel_',int2str(ifile),'.mat'];
    end
    
elseif ~isfield(OPTIONS,'rooot') & ~isempty(OPTIONS.HeadModelFile)
    OPTIONS.rooot = strrep(OPTIONS.HeadModelFile,'headmodel.mat','');    
end
if  ~isfield(OPTIONS,'rooot') % Last Chance
    OPTIONS.rooot = 'bst_';
end

if ~isempty(findstr('sphere',[OPTIONS.Method{:}])) % Single or nested-sphere approaches 
    
    [Param([MEGndx, EEGndx, EEGREFndx]).Center]       = deal(OPTIONS.HeadCenter);
    [Param([MEGndx, EEGndx, EEGREFndx]).Radii]        = deal(OPTIONS.Radii);
    [Param([MEGndx, EEGndx, EEGREFndx]).Conductivity] = deal(OPTIONS.Conductivity);
    
    if EEG & strcmpi('eeg_3sphereberg',lower(OPTIONS.Method{DataType.EEG})) % BERG APPROACH
        
        if ~isempty(OPTIONS.HeadModelFile) & exist(OPTIONS.HeadModelFile,'file')
            if OPTIONS.Verbose, bst_message_window('Checking for previous EEG "BERG" parameters'), end
            ParamOld = load(OPTIONS.HeadModelFile,'Param');
            iFlag = 0; % Flag
            if isfield(ParamOld.Param,'Berg') & ~isempty(ParamOld.Param(1).Radii)
                % Check if these older parameters were computed with the same as current radii and conductivity values
                if ParamOld.Param(1).Radii ~= Param(1).Radii;
                    iFlag = 1;
                end
                if ParamOld.Param(1).Conductivity ~= Param(1).Conductivity;
                    iFlag = 1;
                end
            else
                iFlag = 1;
            end
        else
            iFlag = 1;
        end
        
        if iFlag == 1
            if OPTIONS.Verbose , bst_message_window('Computing EEG "BERG" Parameters. . .'), end
            [mu_berg_tmp,lam_berg_tmp] = berg(OPTIONS.Radii,OPTIONS.Conductivity);
            Param(EEGndx(1)).Berg.mu = mu_berg_tmp; clear mu_berg_tmp
            Param(EEGndx(1)).Berg.lam = lam_berg_tmp; clear lam_berg_tmp
            if OPTIONS.Verbose, bst_message_window({'Computing EEG "BERG" Parameters -> DONE',' '}), end
        else
            if OPTIONS.Verbose , bst_message_window('Using Previous EEG "BERG" Parameters'), end
            Param(EEGndx(1)).Berg.mu = ParamOld.Param(1).Berg.mu;
            Param(EEGndx(1)).Berg.lam = ParamOld.Param(1).Berg.lam; clear ParamOld
        end
        [Param.Berg]= deal(Param(EEGndx(1)).Berg);        
        
    end    
    
end

if ~isempty(findstr('bem',[OPTIONS.Method{:}])) % BEM approaches - Compute transfer matrices

    if OPTIONS.BEM.Interpolative==0 & OPTIONS.VolumeSourceGrid %& isempty(OPTIONS.Cortex) % User wants volumic grid : force Interpolative approach
        OPTIONS.BEM.Interpolative = 1;    % CBB (SB, 07-May-2004)| Should work also for volumic grid
        hwarn = warndlg('Volumic Source Grid BEM is only available for interpolative BEM. BEM computation will be now forced to interpolative. If you want a non-interpolative BEM on a cortical image grid, first uncheck the Volumic Grid box from headmodeler gui.','Limitation from current BrainStorm version');
        drawnow
        waitfor(hwarn)
    end
    
    if MEG
        if ~isempty(findstr('bem',OPTIONS.Method{DataType.MEG})) % BEM is requested for MEG
            BEMChanNdx{DataType.MEG} = MEGndx;
        end
    end
    
    if EEG 
        if ~isempty(findstr('bem',OPTIONS.Method{DataType.EEG})) % BEM is requested for EEG
            BEMChanNdx{DataType.EEG} = sort([EEGndx,EEGREFndx]); % EEGREFndx = [] is average ref
        end
    end
    
    OPTIONS.Param = Param;
    
    if OPTIONS.BEM.Interpolative % Computation of gain matrix over 3D interpolative grid 
        if OPTIONS.BEM.ForceXferComputation
            BEMGaingridFileName = bem_GainGrid(DataType, OPTIONS, BEMChanNdx); % Computation of the BEM gain matrix on the 3D interpolative grid for MEG and/or EEG data
        else
            try 
                BEMGaingridFileName = OPTIONS.BEM.GaingridFileName;
            catch
                cd(User.STUDIES)
                BEMGaingridFileName = bem_GainGrid(DataType, OPTIONS, BEMChanNdx); % Computation of the BEM gain matrix on the 3D interpolative grid for MEG and/or EEG data
            end
        end
   
        if MEG & EEG
            [Param(MEGndx).bem_gaingrid_mfname] = deal(BEMGaingridFileName.MEG);
            [Param([EEGndx]).bem_gaingrid_mfname] = deal(BEMGaingridFileName.EEG);
        else
            [Param(:).bem_gaingrid_mfname] = deal(BEMGaingridFileName);
        end    
        
    else % BEM gain matrix computation over cortical surface
        
        
        BEMGaingridFileName = bem_GainGrid(DataType, OPTIONS, BEMChanNdx); 
        
        [Param(:).bem_gaingrid_mfname] = deal('');
        
    end

  
    OPTIONS = rmfield(OPTIONS,'Param');
    
    ndx = sort([BEMChanNdx{:}]);
    
    
    [Param(ndx).Center] = deal([]);
    test=0; % this test is nowhere defined ?
    
    if OPTIONS.VolumeSourceGrid
        % Now define the outer scalp envelope for the source volume gridding if requested
          
        if OPTIONS.BEM.ForceXferComputation | ~test
            global nfv
            SubjectFV.vertices = nfv(end).vertices'; 
            SubjectFV.faces = nfv(end).faces; clear Faces
        else
            load(fullfile(User.SUBJECTS,fileparts(OPTIONS.Subject),OPTIONS.BEM.EnvelopeNames{end}.TessFile))
        
            idScalp = find(strcmpi(OPTIONS.BEM.EnvelopeNames{end}.TessName,Comment)); clear Comment
            if isempty(idScalp)
                errodlg(sprintf(...
                    'Scalp tessellation %s was not found in %s.',OPTIONS.BEM.EnvelopeNames{end}.TessName, OPTIONS.BEM.EnvelopeNames{end}.TessFile),...
                    'Error during BEM computation')
                return
            end
            
            SubjectFV.vertices = Vertices{idScalp}'; clear Vertices
            SubjectFV.faces = Faces{idScalp}; clear Faces
        end
        
    end
    
end % if BEM

if EEG
    
    switch(lower(OPTIONS.Method{DataType.EEG}))
    case 'eeg_sphere'
        [Param(EEGndx).EEGType] = deal('EEG_SINGLE');
    case 'eeg_3sphere'
        [Param(EEGndx).EEGType] = deal('EEG_3SHELL');
    case 'eeg_3sphereberg'
        [Param(EEGndx).EEGType] = deal('EEG_BERG');
    case 'eeg_os'
        [Param(EEGndx).EEGType] = deal('EEG_OS');
    case 'eeg_bem'
        [Param(EEGndx).EEGType] = deal('BEM');
        if EEG & ~MEG
            [Param(EEGndx).bem_gaingrid_mfname] = deal(BEMGaingridFileName);    
        end
    end
    
else
    
    tmp = [];
    [Param.Conductivity] = deal(tmp);
    [Param.Berg] = deal(tmp);
    [Param.EEGType] = deal(tmp);   
    
end


%_________________________________________________________________________________________________________________________________________________________________________
% 
%  The following part is an adaptation of the former HEADMODEL_MAKE script  
%
% ________________________________________________________________________________________________________________________________________________________________________


% Allocate function names depending on the selected forward approaches specified in the Method argument
Function = cell(1,length(Channel)); %allocate function names

if MEG
    if ~isempty(findstr('bem',OPTIONS.Method{DataType.MEG})) % MEG BEM
        Function(MEGndx) = deal({'meg_bem'});
        if isfield(Channel(MEGndx(1)),'irefsens')
            Function(Channel(MEGndx(1)).irefsens) = deal({'meg_bem'});
        end
        
    else
        Function(MEGndx) = deal({'os_meg'});
        if isfield(Channel(MEGndx(1)),'irefsens')
            Function(Channel(MEGndx(1)).irefsens) = deal({'os_meg'});
        end
        if isfield(OPTIONS,'BEM')
            OPTIONS = rmfield(OPTIONS,'BEM');
        end
        
    end
end


if EEG
    if ~isempty(findstr('bem',OPTIONS.Method{DataType.EEG})) % MEG BEM
        Function(EEGndx) = deal({'eeg_bem'});
    else
        Function(EEGndx) = deal({'eeg_sph'});
        if isfield(OPTIONS,'BEM')
            OPTIONS = rmfield(OPTIONS,'BEM');
        end
    end
    
end            
% --------------------------------------------------------------------------------------------------------

DIMS = [3 12]; % number of columns for each parametric source model: Current Dipole / Current Multipole
HeadModel.Param = Param; 
HeadModel.Function = Function;

if OPTIONS.VolumeSourceGrid % if SearchGain matrices are requested
    
    SearchGain = cell(1,6); % One cell for each source order + pairs of sync. sources - keep 2 'phantom' cells for former Magnetic Dipole model
    
    %------------------------------ Computing SearchGridLoc and SearchGain-------------------------------------------
    
    if ~isempty(findstr(MegMethod,'bem')) | (EEG & ~MEG)
        if OPTIONS.Verbose, bst_message_window(...
                {'Forward modeling not available for EEG''s Current Mulipole models or the BEM approach','Computing forward fields for Current Dipole models. . .'...
                    ,' '}), end
        %return
        OPTIONS.SourceModel = [-1]; % Force current dipole source model when MEG BEM or any EEG is requested
    end
    
    if OPTIONS.Verbose, bst_message_window(...
            'Computing the Gain Matrices from the Volumic Search Grid. . .'), end
    
    % Prepare call to source_grids
    
    GUI.VALIDORDER = OPTIONS.SourceModel;
    GUI.SPACING = 10*OPTIONS.VolumeSourceGridSpacing; % Pass it to source_grids in millimeters
    GUI.VERBOSE = OPTIONS.Verbose;
    GUI.MEG = MEG;
    GUI.EEG = EEG;
    
    if ~exist('SubjectFV','var')
        % Create a spherical tessellation of the INNER SKULL surface boundary for source gridding 
        
        [x,y,z] = sphere(30);
        %Scale to Radius BrainStorm.R
        [TH,PHI,R] = cart2sph(x,y,z); 
        if OPTIONS.Radii(1) == 0 
            error('The sphere radius you have specified is ZERO - please change it to a non null positive value')
        end
        R = OPTIONS.Radii(1) * ones(size(R));
        [x,y,z] = sph2cart(TH,PHI,R);
        
        % Tessellate both hemispheres     
        Im = find(z>=0);
        tri = delaunay(x(Im),y(Im));
        
        %Gather everything into the same patch
        SubjectFV.vertices = [x(Im),y(Im),z(Im);...
                x(Im),y(Im),-z(Im)];
        
        % Translate about the head center
        SubjectFV.vertices = SubjectFV.vertices + repmat(OPTIONS.HeadCenter',size(SubjectFV.vertices,1),1);
        
        SubjectFV.faces = [tri;tri(:,[3 2 1])+max(tri(:))];
        
        clear x y z R TH PHI Im tri
        
    end
    
    HeadModel = source_grids(HeadModel,Channel,SubjectFV,GUI); 
    clear SubjectFV
    
    if OPTIONS.Verbose, bst_message_window({...
                '-> DONE',...
                'Now saving HeadModel file. . .'...
            })
    end
    
    % Now save VolumeSourceGrid headmodel(s)
    SaveHeadModel.Param = HeadModel.Param;
    SaveHeadModel.Function = HeadModel.Function;
    
    if strcmpi(OPTIONS.HeadModelName,'Default') % Specify default HeadModelName 
        if MEG & EEG
            SaveHeadModel.HeadModelName = sprintf('%s | %s', OPTIONS.Method{DataType.MEG},OPTIONS.Method{DataType.EEG});
        elseif MEG
            SaveHeadModel.HeadModelName = sprintf('%s', OPTIONS.Method{DataType.MEG});
        elseif EEG
            SaveHeadModel.HeadModelName = sprintf('%s', OPTIONS.Method{DataType.EEG});
        end
    else
        SaveHeadModel.HeadModelName = OPTIONS.HeadModelName;
    end
    
    if ~isempty(OPTIONS.HeadModelFile), 
        [HeadModelFilePath,HeadModelFile,ext] = fileparts(OPTIONS.HeadModelFile);
    end
    
try
    cd(User.STUDIES)
catch
end

    
    for k = 1:length(GUI.VALIDORDER) % One file per source order 
        
        % Specify Source Model Name and Other Fields 
        if GUI.VALIDORDER(k) == -1
            SourceOrderString = 'CD'; % Current Dipole
        elseif GUI.VALIDORDER(k) == 0
            SourceOrderString = 'MD'; % Magnetic Dipole
        elseif GUI.VALIDORDER(k) == 1
            SourceOrderString = 'CME'; % Current Multipole
        end
        
        [SaveHeadModel.Param.Order] = deal(GUI.VALIDORDER(k));
        SaveHeadModel.SourceOrder = GUI.VALIDORDER(k); % Alternative to previous line
        SaveHeadModel.HeadModelType = 'SearchGrid';
        if MEG 
            SaveHeadModel.MEGMethod = OPTIONS.Method{DataType.MEG};
        end
        if EEG 
            SaveHeadModel.EEGMethod = OPTIONS.Method{DataType.EEG};
        end
        
        SaveHeadModel.GridName = {sprintf('%s Volumic Grid : %3.1f cm spacing',SourceOrderString, OPTIONS.VolumeSourceGridSpacing)};
        SaveHeadModel.GridLoc  = HeadModel.GridLoc(GUI.VALIDORDER(k)+2);
        SaveHeadModel.Gain = {single(HeadModel.Gain{GUI.VALIDORDER(k)+2})}; % Convert to single precision gain matrix as in BST MMII convention
        
        % now collect together and save
        if ~isempty(OPTIONS.HeadModelFile), % User has provided a headmodel file name 
            SaveHeadModelFile = fullfile(HeadModelFilePath,[HeadModelFile,...
                    sprintf('VolGrid_%s',SourceOrderString),...
                    ext]);
            ifile = 0;
            while exist(SaveHeadModelFile,'file') % Do not overwrite headmodel files
                ifile = ifile + 1;
                SaveHeadModelFile = fullfile(HeadModelFilePath,[HeadModelFile,...
                        sprintf('VolGrid_%s',SourceOrderString),...
                        '_',int2str(ifile),ext]);
            end
            
            
            if OPTIONS.Verbose, bst_message_window({...
                    sprintf('Writing HeadModel file:'),...
                    sprintf('%s', SaveHeadModelFile)...
                })
            end
            
            if strcmp(lower(OPTIONS.HeadModelFileOld),'default')
                OPTIONS.HeadModelFile = SaveHeadModelFile;
            end
            
            save_fieldnames(SaveHeadModel, SaveHeadModelFile);
            
            if OPTIONS.Verbose, bst_message_window(...
                    {'-> DONE',' '}), end
        end
        
    end
    
    OPTIONS.HeadModelName = SaveHeadModel.HeadModelName;
    
    if OPTIONS.Verbose, bst_message_window(...
            {'Grid gain matrices are saved','RAP-MUSIC & Least-Squares Fits Approaches are now available',' '})
    end
    
end


%% -------------------------------------------------------------------------
%
%  Now proceed to forward modeling of cortical grids or at some other specific source locations
%
%% -------------------------------------------------------------------------


if ~isempty(OPTIONS.Cortex), % subject has cortical vertices as source supports
    % First make room in memory
    HeadModel.SearchGain = [];
    clear SearchGridLoc SearchGain G
    
    if OPTIONS.Verbose, bst_message_window({...
                'Computing the Image Gain Matrices (this may take a while). . .'})
    end
    
    % Find the cortical grid where to compute the forward model in the tessellation file
    try 
        ImageGrid = load(fullfile(User.SUBJECTS,OPTIONS.Cortex.FileName),'Comment','Vertices');
    catch
        ImageGrid = load(OPTIONS.Cortex.FileName,'Comment','Vertices');
    end
    
    if ~isfield(OPTIONS.Cortex,'iGrid')
        OPTIONS.Cortex.iGrid = 1;
    end
    
    if OPTIONS.Cortex.iGrid > length(ImageGrid.Comment) % Corresponding cortical surface not found
        errordlg(sprintf('Cortical Tessellation file "%s" does not contain %d surfaces',...
            OPTIONS.Cortex.FileName,OPTIONS.Cortex.iGrid))
        return
    end
    
    GridLoc = ImageGrid.Vertices{OPTIONS.Cortex.iGrid}; % keep only the desired cell
    ImageGrid = rmfield(ImageGrid,'Vertices');
    
elseif ~isempty(OPTIONS.SourceLoc) % Specific source locations are provided
    if size(OPTIONS.SourceLoc,1) ~= 3
        OPTIONS.SourceLoc = OPTIONS.SourceLoc';
    end
    GridLoc = OPTIONS.SourceLoc;
    
else
    rmfield(OPTIONS,'rooot');    
    if nargout == 1
        varargout{1} = OPTIONS;
    else
        varargout{1} = HeadModel;
        varargout{2} = OPTIONS;
    end
    return % No ImageGrid requested
    
end

if ~isempty(OPTIONS.Cortex)
    GridName{OPTIONS.Cortex.iGrid} = ImageGrid.Comment{OPTIONS.Cortex.iGrid}; 
    OPTIONS.Cortex.Name  = ImageGrid.Comment{OPTIONS.Cortex.iGrid};
end

%for i = 1 % CHEAT - Dipoles only here
for Order = OPTIONS.SourceModel % Compute gain matrices for each requested source models    (-1 0 1)
    
    i = 1; % Index to cell in headmodel cell arrays (MMII convention)
    
    switch(Order)
    case -1
        SourceOrderString = 'CD'; % Current Dipole
        Dims = DIMS(1);% number of columns per source
    case 0
        errordlg(sprintf('Unauthorized Source Model Order %d.',iSrcModel),'Wrong HeadModel parameter assignment')
        return
        %SourceOrderString = 'MD'; % Magnetic Dipole  - OBSOLETE
    case 1
        SourceOrderString = 'CME'; % Current Multipole
        Dims = DIMS(2);% number of columns per source
    otherwise
        errordlg(sprintf('Unauthorized Source Model Order %d.',iSrcModel),'Wrong HeadModel parameter assignment')
        return
    end
    
    if OPTIONS.Verbose, 
        bst_message_window(...
            'wrap',...
            sprintf('Computing Gain Matrix for the %s Model',SourceOrderString))
        h = waitbar(0,sprintf('Computing Gain Matrix for the %s Model',SourceOrderString));
        pos = get(h,'position');
    end
    
    if ~isempty(OPTIONS.Cortex) & OPTIONS.ApplyGridOrient % Use cortical grid
        try
            load(OPTIONS.Cortex.FileName,'Faces');
        catch
            Users = get_user_directory;
            load(fullfile(Users.SUBJECTS,OPTIONS.Cortex.FileName),'Faces');
        end
        
        Faces = Faces{OPTIONS.Cortex.iGrid};
        
        ptch = patch('Vertices',GridLoc','Faces',Faces,'Visible','off');
        set(get(ptch,'Parent'),'Visible','off')
        clear Faces
        GridOrient{i} = get(ptch,'VertexNormals')'; % == {1} as of MMII conventions. Most data in HeadModel files are single-cell cell arrays.
        % Cell array structure kept for backward compatibility with BsT2000.
        delete(ptch);
    
    end
    
    if OPTIONS.ApplyGridOrient % Take cortex normals into account
        if isempty(OPTIONS.GridOrient) & isempty(OPTIONS.SourceOrient) % Consider the cortical patch's normals
            
            if isempty(OPTIONS.Cortex)  % Specific source locations in .SourceLoc but nohing in. SourceOrient
                OPTIONS.ApplyGridOrient = 0; % No source orientation specified: Force computation of full gain matrix         
            else
                [nrm,GridOrient{i}] = colnorm(GridOrient{i});
                % Now because some orientations may be ill-set to [0 0 0] with the get(*,'VertexNormals' command) 
                % set these orientation to arbitrary [1 1 1]:
                izero = find(nrm == 0);clear nrm
                if ~isempty(izero)
                    GridOrient{i}(:,izero) = repmat([1 1 1]'/norm([1 1 1]),1,length(izero));
                end
                clear izero
                
            end
            
        elseif ~isempty(OPTIONS.GridOrient) % Apply user-defined cortical source orientations
            
            if size(OPTIONS.GridOrient,2) == size(GridLoc,2) % Check size integrity
                GridOrient{i} = OPTIONS.GridOrient;
                [nrm,GridOrient{i}] = colnorm(GridOrient{i});
                clear nrm
            else
                errordlg(sprintf('The source orientations you have provided are for %0.f sources. Considered cortical surface has %0.f sources. Computation aborted',...
                    size(OPTIONS.GridOrient,2),size(GridOrient{i},2)));
                return
            end
            
        elseif ~isempty(OPTIONS.SourceOrient) % Apply user-defined specific source orientations
            
            if size(OPTIONS.SourceOrient,2) == size(GridLoc,2) % Check size integrity
                GridOrient{i} = OPTIONS.SourceOrient;
                [nrm,GridOrient{i}] = colnorm(GridOrient{i});
                clear nrm
            else
                errordlg(sprintf('The source orientations you have provided are for %0.f sources. Computation aborted',...
                    size(OPTIONS.SourceOrient,2)));
                return
            end
            
        end
    end
    
    nv = size(GridLoc,2); % number of grid points
    
    if ~isempty(OPTIONS.ImageGridFile) % Save cortical gain matrix in a binary file
        
        [PATH,NAME,EXT,VER] = fileparts(OPTIONS.HeadModelFile);
        
        if strcmp(lower(OPTIONS.ImageGridFile),'default') 
            if exist('GridName','var') % Cortical support was specified
                destname = [NAME,'_Gain_',strrep(ImageGrid.Comment{OPTIONS.Cortex.iGrid},' ',''),'_',SourceOrderString,'.bin']; % New naming (March, 19 - 2002)
                k = 1;
                try
                    while exist(destname,'file') % Don't write over existing .bin gain matrix file
                        destname = [NAME,'_Gain_',strrep(ImageGrid.Comment{OPTIONS.Cortex.iGrid},' ',''),'_',SourceOrderString,'_',int2str(k),'.bin']; % New naming (March, 19 - 2002)
                        k = k+1;
                    end
                catch
                    while exist(destname,'file') % Don't write over existing .bin gain matrix file
                        destname = [NAME,'_Gain_',strrep(ImageGrid.Comment{OPTIONS.Cortex.iGrid},' ',''),'_',SourceOrderString,'_',int2str(k),'.bin']; % New naming (March, 19 - 2002)
                        k = k+1;
                    end
                end
                
                clear k    
                
            else  % Specific location was provided in .SourceLoc
                destname = [NAME,'_Gain_SpecLoc_',SourceOrderString,'.bin']; % Specific location was provided in .SourceLoc
                k = 1;
                while exist(destname,'file') % Don't write over existing .bin gain matrix file
                    destname = [NAME,'_Gain_SpecLoc_',SourceOrderString,'_',int2str(k),'.bin']; 
                    k = k+1;
                end
                clear k
            end
            OPTIONS.ImageGridFile = destname;
        else
            [PATH,destname,ext,ver] = fileparts(OPTIONS.ImageGridFile);
            destname = [destname, ext];
        end
        
        % Check whether this name exists - if yes, don't overwrite
        try
            cd(fullfile(User.STUDIES,PATH))
        catch
            cd(PATH)
        end
        
        hdml = 1;
        while exist(destname,'file')
            if hdml == 1
                destname = strrep(destname,'.bin',['_',int2str(hdml),'.bin']);
            else
                destname = strrep(destname,[int2str(hdml-1),'.bin'],[int2str(hdml),'.bin']);
            end
            
            hdml = hdml+1;
        end
        clear hdml
        
        destnamexyz = [destname(1:end-4) '_xyz.bin'];
        fdest = fopen(destname,'w','ieee-be');
        fdestxyz = fopen(destnamexyz,'w','ieee-be');
        if ((fdest < 0) | (fdestxyz < 0))
            errordlg('Error creating the file for the cortical image forward model; Please check for disk space availability')
            return
        end
        frewind(fdest);
        frewind(fdestxyz)
        
        fwrite(fdest,length([OPTIONS.Channel]),'uint32'); 
        fwrite(fdestxyz,length([OPTIONS.Channel]),'uint32'); 
        
        %         % BEM and non-interpolative, store gain matrix already computed 
        %         if isfield(OPTIONS,'BEM') % BEM computation
        %             if ~OPTIONS.BEM.Interpolative
        %                 global GBEM_grid
        %             end
        %         end
        
        if ((fdest < 0) | (fdestxyz < 0)) , errordlg('Please Check Write Permissions', ['Cannot create ',destname]), return, end
    else % If no ImageGridFile was specified
        %OPTIONS.ImageGridBlockSize = nv; % Force a one-time computation of all source forward fields
    end
    
    src_ind = 0;
    jj = 0; % Number of OPTIONS.ImageGridBlockSize 
    
    
    for j = 1:(OPTIONS.ImageGridBlockSize):nv,  
        jj = jj+1;
        if 1%OPTIONS.ApplyGridOrient
            ndx = [0:OPTIONS.ImageGridBlockSize-1]+j;
        else
            ndx = [0:3*OPTIONS.ImageGridBlockSize-1]+j;
        end
        
        if 1%OPTIONS.ApplyGridOrient
            if(ndx(end) > nv),  % last OPTIONS.ImageGridOPTIONS.ImageGridBlockSizeSize too long
                ndx = [ndx(1):nv];
            end
        else
            if(ndx(end) > 3*nv),  
                ndx = [ndx(1):3*nv];
            end
        end
        
        % Compute MEG
        if MEG & ~isempty(MEGndx) 
            Gmeg = NaN*zeros(length(MEGndx),Dims*length(ndx)); 
            if ~isempty(Function{MEGndx(1)})
                
                if MEG & EEG
                    clear('gain_bem_interp2'); % Free persistent variables to avoid confusion
                end
                
                if isfield(OPTIONS,'BEM') % BEM computation
                    if ~OPTIONS.BEM.Interpolative % ~interpolative : retrieve stored gain matrix 
                        global GBEM_grid
                        tmpndx = [3*(ndx(1)-1)+1:min([3*nv,3*ndx(end)])];%[0:3*OPTIONS.ImageGridBlockSize-1]+j;
                        Gmeg = GBEM_grid(MEGndx,tmpndx);
                    end
                else
                    Gmeg = feval(Function{MEGndx(1)},GridLoc(:,ndx),Channel,Param,Order,OPTIONS.Verbose);
                end

            end
        else
            Gmeg = []; 
        end
        
        if EEG & ~isempty(EEGndx) & i==1 % % Order -1 only for now in EEG 
            
            Geeg = NaN*zeros(length(EEGndx),Dims*length(ndx)); 
            if ~isempty(Function{EEGndx(1)})
                if MEG & EEG
                    clear('gain_bem_interp2'); % Free persistent variables to avoid confusion
                end

                if isfield(OPTIONS,'BEM') % BEM computation
                    if ~OPTIONS.BEM.Interpolative % ~interpolative : retrieve stored gain matrix 
                        global GBEM_grid
                        %Geeg = GBEM_grid(EEGndx,[j:min([3*nv,j+3*OPTIONS.ImageGridBlockSize-1])]);
                        tmpndx = [3*(ndx(1)-1)+1:min([3*nv,3*ndx(end)])];%[0:3*OPTIONS.ImageGridBlockSize-1]+j;
                        %min(tmpndx ), max(tmpndx)
                        Geeg = GBEM_grid(EEGndx,tmpndx);
                    end
                else
                    Geeg = feval(Function{EEGndx(1)},GridLoc(:,ndx),Channel,Param,Order,OPTIONS.Verbose);
                end
                
            end
        else
            Geeg = [];
        end
        
            G = NaN*zeros(length(OPTIONS.Channel),length(ndx));
            Gxyz = NaN*zeros(length(OPTIONS.Channel),Dims*length(ndx));
        
        if OPTIONS.Verbose, bst_message_window(...
                ['Computing Cortical Gain Vectors. . . Block # ',int2str(jj),...
                    ' of ',int2str(length(1:OPTIONS.ImageGridBlockSize:nv))])
            
            hh = waitbar(0,['Computing Cortical Gain Vectors. . . Block # ',int2str(jj),...
                    ' of ',int2str(length(1:OPTIONS.ImageGridBlockSize:nv))]);
            set(hh,'Position',[pos(1), pos(2)+pos(4),pos(3),pos(4)])
            drawnow
        end
        
        src = 0;
        % Options on cortical source orientation
%        if OPTIONS.ApplyGridOrient % Apply source orientation
            
            for k = 1:Dims:Dims*length(ndx)-2
                src = src+1;
                src_ind = src_ind+1;
                if ~isempty(Gmeg)
                    G(MEGndx,src) = Gmeg(:,k:k+2) * GridOrient{1}(:,src_ind);   
                end
                
                if ~isempty(Geeg)&(i==1)% % Order -1 only 
                    G(EEGndx,src) = Geeg(:,k:k+2) * GridOrient{1}(:,src_ind);   
                end
                if OPTIONS.Verbose, 
                    if ~rem(src,5000)
                        waitbar(src/length(1:Dims:Dims*length(ndx)-2),hh)
                    end
                end
                
            end
            
%        else % Do not apply cortical orientation. Keep full gain matrix at each cortical location
            
            if MEG 
                Gxyz(MEGndx,:) = Gmeg;   
            end
            if EEG
                Gxyz(EEGndx,:) = Geeg;
            end
            
%        end
        
        if OPTIONS.Verbose, 
           if ~rem(src_ind,1000)
               waitbar(src_ind/nv,h)
           end
           
            delete(hh)
        end
        
        clear Gmeg Geeg
        
        if ~isempty(OPTIONS.ImageGridFile) 
            % 4 bytes per element, find starting point
            offset = 4 + (ndx(1)-1)*length(Channel)*4;
            
            % Matlab 6.5.0 bug, Technical Solution Number: 1-19NOK
            % http://www.mathworks.com/support/solutions/data/1-19NOK.html?solution=1-19NOK
            % must REWIND first before fseek.
            % Apparently fixed in 6.5.1
            % JCM 27-May-2004
            status = fseek(fdest,0,'bof');
            status = fseek(fdest,offset,'bof');
            if(status == -1),
                errordlg('Error writing Image Gain Matrix file'); return
            end

            fwrite(fdest,G,'float32');

            %save xyz forward matrix
            offset = 4 + 3*(ndx(1)-1)*length(Channel)*4;
            status = fseek(fdestxyz,0,'bof');
            status = fseek(fdestxyz,offset,'bof');
            if(status == -1),
                errordlg('Error writing Image Gain Matrix file'); return
            end
            fwrite(fdestxyz,Gxyz,'float32');
        
        end
        
    end
    if exist('destname','var') & ~isempty(OPTIONS.Cortex)
        Gain{OPTIONS.Cortex.iGrid}{i} = destname;
        OPTIONS.ImageGridFile = destname;    
    end
    
    if OPTIONS.Verbose
        close(h)
    end
    
    if ~isempty(OPTIONS.ImageGridFile) 
        fclose(fdest);
        if OPTIONS.Verbose, bst_message_window({...
                    sprintf('Computing Gain Matrix for the %s Model -> DONE',SourceOrderString),...
                    sprintf('Saved in:'),...
                    sprintf('%s',destname)...
                    }),
        end
    end
    
    
    % Now save ImageGrid headmodel(s)-----------------------------------------------------------------------------------
    OPTIONS.HeadModelFile = OPTIONS.HeadModelFileOld; % Use original HeadModelFile entry
    if strcmp(lower(OPTIONS.HeadModelFile),'default')
        if ~isfield(OPTIONS,'rooot')
            OPTIONS.rooot = OPTIONS.PrefixFileName;
        end
        OPTIONS.HeadModelFile = [OPTIONS.rooot,'headmodel.mat'];
        ifile = 0; 
        while exist(OPTIONS.HeadModelFile,'file') % Do not overwrite headmodel files
            ifile = ifile + 1;
            OPTIONS.HeadModelFile = [OPTIONS.rooot,'headmodel_',int2str(ifile),'.mat'];
        end
        
    elseif ~isfield(OPTIONS,'rooot') & ~isempty(OPTIONS.HeadModelFile)
        OPTIONS.rooot = strrep(OPTIONS.HeadModelFile,'headmodel.mat','');    
    end
    
    SaveHeadModel.Param = HeadModel.Param;
    SaveHeadModel.Function = HeadModel.Function;
    
    if MEG 
        SaveHeadModel.MEGMethod = OPTIONS.Method{DataType.MEG};
    end
    if EEG 
        SaveHeadModel.EEGMethod = OPTIONS.Method{DataType.EEG};
    end
    
    if strcmpi(OPTIONS.HeadModelName,'Default') % Specify default HeadModelName 
        if MEG & EEG
            SaveHeadModel.HeadModelName = sprintf('%s | %s', OPTIONS.Method{DataType.MEG},OPTIONS.Method{DataType.EEG});
        elseif MEG
            SaveHeadModel.HeadModelName = sprintf('%s', OPTIONS.Method{DataType.MEG});
        elseif EEG
            SaveHeadModel.HeadModelName = sprintf('%s', OPTIONS.Method{DataType.EEG});
        end
    else
        SaveHeadModel.HeadModelName = OPTIONS.HeadModelName;    
    end
    
    if ~isempty(OPTIONS.HeadModelFile), % and is not empty
        [HeadModelFilePath,HeadModelFile,ext] = fileparts(OPTIONS.HeadModelFile);
    end
    
    % Assign proper GridName
    switch(Order)
    case -1
        SourceOrderString = 'CD'; % Current Dipole
    case 0
        SourceOrderString = 'MD'; % Magnetic Dipole
    case 1
        SourceOrderString = 'CME'; % Current Multipole
    end
    
    [SaveHeadModel.Param.Order] = deal(Order);
    SaveHeadModel.SourceOrder = Order; % Alternative to previous line
    SaveHeadModel.HeadModelType = 'ImageGrid';
    
    
    if ~isempty(OPTIONS.Cortex)
        SaveHeadModel.GridName = {sprintf('%s Surface Grid : %s',SourceOrderString, OPTIONS.Cortex.Name)};
        if ~isempty(OPTIONS.ImageGridFile) 
            SaveHeadModel.Gain = {OPTIONS.ImageGridFile}; % Store Image Grid File name in the headmodel.mat file
        end
        
        SaveHeadModel.GainCovar{1} = cell(1); % May compute it later within headmodel_make
        SaveHeadModel.GainCovarName ='';
        SaveHeadModel.GridLoc{1} = OPTIONS.Cortex.FileName;
        if 1%OPTIONS.Cortex.iGrid > 1 % Meaning that cortical support is not the fisrt surface in tessellation file (default in MMII)
            % Add yet another field to headmodel file to specify this information.
            SaveHeadModel.iGrid = OPTIONS.Cortex.iGrid;
        end
    else
        SaveHeadModel.GainCovar = [];
        SaveHeadModel.GainCovarName = [];
        SaveHeadModel.GridName = [];
        SaveHeadModel.GridLoc = {OPTIONS.SourceLoc};
        SaveHeadModel.Gain = {OPTIONS.ImageGridFile};
    end
    if ~OPTIONS.ApplyGridOrient
        SaveHeadModel.GridOrient = [];
    else
        SaveHeadModel.GridOrient = {OPTIONS.SourceOrient};
    end
    
    
    % now collect together and save
    if ~isempty(OPTIONS.HeadModelFile) & exist('destname','var')
        % HeadModel file name
        SaveHeadModelFile = fullfile(HeadModelFilePath,[HeadModelFile,...
                sprintf('SurfGrid_%s',SourceOrderString),...
                ext]);
        ifile = 0;
        while exist(SaveHeadModelFile,'file') % Do not overwrite headmodel files
            ifile = ifile + 1;
            SaveHeadModelFile = fullfile(HeadModelFilePath,[HeadModelFile,...
                    sprintf('SurfGrid_%s',SourceOrderString),...
                    '_',int2str(ifile),ext]);
        end
        if OPTIONS.Verbose, bst_message_window({...
                    sprintf('Writing Cortical Image Support HeadModel file:'),...
                    sprintf('%s', SaveHeadModelFile)...
                })
        end
        
        if strcmp(lower(OPTIONS.HeadModelFileOld),'default')
            OPTIONS.HeadModelFile = SaveHeadModelFile;
        end
        
        try
            save_fieldnames(SaveHeadModel, SaveHeadModelFile);
        catch
            cd(User.STUDIES)
            save_fieldnames(SaveHeadModel, SaveHeadModelFile);
        end
        
        
        if OPTIONS.Verbose, bst_message_window(...
                {'-> DONE',' '}), end
    end
    
    % Save completed -----------------------------------------------------------------------------------
  
end % Cortical gain matrix for each source model


%% -------------------------------------------------------------------------
if isfield(OPTIONS,'rooot')
    OPTIONS = rmfield(OPTIONS,'rooot');
end

if nargout == 0
    clear G SearchGain
elseif nargout == 1
    varargout{1} = OPTIONS;
else
    if OPTIONS.ImageGridBlockSize < nv
        varargout{1} = OPTIONS.HeadModelFile;
    else
        varargout{1} = G;
    end
    
    varargout{2} = OPTIONS;
end

fclose('all');

if OPTIONS.Verbose, bst_message_window(...
        {'The head model has been properly designed and written to disk'})
end

if OPTIONS.Verbose, bst_message_window({...
            sprintf(...
            'Head modeling ends (%s - %dH %dmin %ds (took %3.2f seconds))',date,clockk(4:6), etime(clock,time0)),...
            '__________________________________________'...
        })
end

%------------------------------------------------------------------------------------------------------------------------------
% 
%                                         SUB-FUNCTIONS
% 
%------------------------------------------------------------------------------------------------------------------------------

function BEMGaingridFname = bem_GainGrid(DataType,OPTIONS,BEMChanNdx)

% Computation of the BEM gain matrix on the 3D interpolative grid for MEG and/or EEG data
% DataType : a structure with fields MEG and EEG. DataType.MEG (res. DataType.EEG) is set to 1 if MEG (res. EEG) data is available
% OPTIONS : the OPTIONS structure, input argument from bst_headmodeler
% BEMChanNdx :  a cell array of channel indices such that : BEMChanNdx{DataType.EEG} = EEGndx (res. MEG).
%
% BEMGaingridFname: a structure with fields:
%                 .EEG : string with the name of the file containing the gain matrix of the 3D BEM interpolative grid in EEG
%                 .MEG : same as .EEG respectively to MEG BEM model.


% Detect the requested BEM computations: MEG and/or EEG___________________________________

User = get_user_directory;
if isempty(User)
    User.STUDIES  = pwd;
    User.SUBJECTS = pwd;
end

MEG = ~isempty(BEMChanNdx(DataType.MEG)); %  == 1 if MEG is requested 
EEG = ~isempty(BEMChanNdx(DataType.EEG)); %  == 1 if EEG is requested
if MEG 
    MEGndx = BEMChanNdx{DataType.MEG};
end
if EEG 
    %EEGndx = BEMChanNdx{DataType.EEG};
    EEGndx = OPTIONS.EEGndx; % EEG sensors (not including EEG reference channel, if any)
    EEGREFndx = good_channel(OPTIONS.Channel,[],'EEG REF');
end

if MEG & EEG
    [Param(:).mode] = deal(3);
elseif ~MEG & EEG
    [Param(:).mode] = deal(1);
elseif MEG & ~EEG
    [Param(:).mode] = deal(2);
else
    errordlg('Please check that the method requested for forward modeling has an authorized name');
    return
end

% BEM parameters ________________________________________________________________________
% Determine what basis functions to use
constant = ~isempty(strmatch('constant',lower(OPTIONS.BEM.Basis)));
if constant == 0
    [Param(:).basis_opt] = deal(1);
else
    [Param(:).basis_opt] = deal(0);
end

% Determine what Test to operate
collocation = ~isempty(strmatch('collocation',lower(OPTIONS.BEM.Test)));
if collocation == 0
    [Param(:).test_opt] = deal(1);
else
    [Param(:).test_opt] = deal(0);
end

% Insulated-skull approach
isa = OPTIONS.BEM.ISA;
if isa == 1
    [Param(:).ISA] = deal(1);
else
    [Param(:).ISA] = deal(0);
end

Param.Ntess_max = OPTIONS.BEM.NVertMax;
Param.Conductivity = deal(OPTIONS.Conductivity);

%_________________________________________________________________________________________

% Load surface envelopes information______________________________________________________

% Find the indices of the enveloppes selected for BEM computation
if ~isfield(OPTIONS.BEM,'EnvelopeNames')
    errordlg('Please specify the ordered set of head-tissue envelopes by filling the OPTIONS.BEM.EnvelopeNames field')
    return
end
if isempty(OPTIONS.BEM.EnvelopeNames)
    errordlg('Please specify the ordered set of head-tissue envelopes by filling the OPTIONS.BEM.EnvelopeNames field')
    return
end
for k = 1:length(OPTIONS.BEM.EnvelopeNames)

    try
        load(fullfile(User.SUBJECTS,OPTIONS.subjectpath,OPTIONS.BEM.EnvelopeNames{k}.TessFile),'Comment');
    catch % Maybe user is using command line call to function with absolute-referenced files OPTIONS.*.TessFile
        try
            OPTIONS.BEM.EnvelopeNames{k}.TessFile = [OPTIONS.BEM.EnvelopeNames{k}.TessFile,'.mat'];
            load(fullfile(User.SUBJECTS,OPTIONS.subjectpath,OPTIONS.BEM.EnvelopeNames{k}.TessFile),'Comment');
        catch
            cd(User.SUBJECTS)
            load(OPTIONS.BEM.EnvelopeNames{k}.TessFile,'Comment');
        end
    end

    Comment = strrep(Comment,' ','');
    % find surface in current tessellation file
    OPTIONS.BEM.EnvelopeNames{k}.SurfId = find(strcmpi(OPTIONS.BEM.EnvelopeNames{k}.TessName,Comment));
    IDs(k) = OPTIONS.BEM.EnvelopeNames{k}.SurfId;
    if isempty(OPTIONS.BEM.EnvelopeNames{k}.SurfId)
        errordlg(...
            sprintf('Surface %s was not found in file %s',...
            OPTIONS.BEM.EnvelopeNames{k}.TessName, OPTIONS.BEM.EnvelopeNames{k}.TessFile))
        return
    end

    % Load vertex locations
    try
        tmp = load(fullfile(User.SUBJECTS,OPTIONS.subjectpath,OPTIONS.BEM.EnvelopeNames{k}.TessFile),'Vertices');
    catch % Maybe user is using command line call to function with absolute-referenced files OPTIONS.*.TessFile
        tmp = load(OPTIONS.BEM.EnvelopeNames{k}.TessFile,'Vertices');
    end
    Vertices{k} = tmp.Vertices{OPTIONS.BEM.EnvelopeNames{k}.SurfId}';

    % Load faces
    try
        tmp = load(fullfile(User.SUBJECTS,OPTIONS.subjectpath,OPTIONS.BEM.EnvelopeNames{k}.TessFile),'Faces');
    catch% Maybe user is using command line call to function with absolute-referenced files OPTIONS.*.TessFile
        tmp = load(OPTIONS.BEM.EnvelopeNames{k}.TessFile,'Faces');
    end

    Faces(k) = tmp.Faces(OPTIONS.BEM.EnvelopeNames{k}.SurfId);
end
clear Comment tmp

cd(User.STUDIES)


%_________________________________________________________________________________________


% Channel Parameters______________________________________________________________________
if MEG
    R_meg1 = zeros(length(MEGndx),3);
    O_meg1 = R_meg1;
    R_meg2 = R_meg1;
    O_meg2 = O_meg1;
    flaggrad = zeros(length(MEGndx),1);% if = 1 - Flag to indicate there are some gradiometers here

    i = 0;
    for k = MEGndx
        i = i+1;
        R_meg1(i,:) = OPTIONS.Channel(k).Loc(:,1)';
        O_meg1(i,:) = OPTIONS.Channel(k).Orient(:,1)';

        if size(OPTIONS.Channel(k).Loc,2) == 2
            if sum(OPTIONS.Channel(k).Loc(:,1)-OPTIONS.Channel(k).Loc(:,2))~=0 % Gradiometer
                R_meg2(i,:) = OPTIONS.Channel(k).Loc(:,2)';
                O_meg2(i,:) = OPTIONS.Channel(k).Orient(:,2)';
                flaggrad(k-min(MEGndx)+1) = 1;
            end
        end

    end
    O_meg1 = (O_meg1' * inorcol(O_meg1'))';
    if exist('O_meg2','var')
        O_meg2 = (O_meg2' * inorcol(O_meg2'))';
    end

    % Handle MEG reference channels if necessary
    %     if isfield(OPTIONS.Channel(MEGndx(1)),'irefsens')
    %         irefsens = OPTIONS.Channel(MEGndx(1)).irefsens;
    %     else
    %         irefsens = [];
    %     end
    irefsens = good_channel(OPTIONS.Channel,[],'MEG REF');

    if ~isempty(irefsens)
        flaggrad_REF = zeros(length(irefsens),1);% if = 1 - Flag to indicate there are some gradiometers here
        R_meg_REF = zeros(length(irefsens),3);
        O_meg_REF = R_meg_REF;
        R_meg_REF = R_meg_REF;
        O_meg_REF = R_meg_REF;

        if ~isempty(irefsens) & ~isempty(OPTIONS.Channel(MEGndx(1)).Comment) % Reference Channels are present
            if OPTIONS.Verbose, bst_message_window('Reference Channels have been detected.'), end
        end
        i = 0;
        for k = irefsens
            i = i+1;
            R_meg_REF(i,:) = OPTIONS.Channel(k).Loc(:,1)';
            O_meg_REF(i,:) = OPTIONS.Channel(k).Orient(:,1)';

            if size(OPTIONS.Channel(k).Loc,2) == 2
                if sum(OPTIONS.Channel(k).Loc(:,1)-OPTIONS.Channel(k).Loc(:,2))~=0 % Reference Gradiometer
                    R_meg_REF2(i,:) = OPTIONS.Channel(k).Loc(:,2)';
                    O_meg_REF2(i,:) = OPTIONS.Channel(k).Orient(:,2)';
                    flaggrad_REF(k-min(irefsens)+1) = 1;
                end
            end

        end

    else
        R_meg_REF = [];
        O_meg_REF = R_meg_REF;
        R_meg_REF = R_meg_REF;
        O_meg_REF = R_meg_REF;
    end

    MEGndx_orig = MEGndx;
else
    R_meg1 = zeros(length(EEGndx),3); % Use dummy channel locations
    O_meg1 = R_meg1;
    R_meg2 = R_meg1;
    O_meg2 = O_meg1;
end

if EEG
    R_eeg = zeros(length([EEGREFndx,EEGndx]),3);
    i = 0;
    for k = [EEGREFndx,EEGndx]
        i = i+1;
        R_eeg(i,:) = OPTIONS.Channel(k).Loc(:,1)';
    end
    if ~MEG
        flaggrad = [];
        irefsens = [];
    end

else
    R_eeg = NaN * zeros(size(R_meg1)); % Dummy coordinates
end

%_________________________________________________________________________________________


%Compute Transfer Matrices__________________________________________________________________

BrainStorm.iscalp = IDs(end); % Index of the outermost surface (scalp, supposedly)

eeg_answ = '';
meg_answ = '';

% Check whether some transfer-matrix files already exist for current study

fn_eeg = sprintf('%s_eegxfer_%s_%s.mat',OPTIONS.rooot, OPTIONS.BEM.Basis,OPTIONS.BEM.Test);
fn_meg = sprintf('%s_megxfer_%s_%s.mat',OPTIONS.rooot, OPTIONS.BEM.Basis,OPTIONS.BEM.Test);
if MEG
    test = exist(fn_meg,'file');
elseif EEG
    test = exist(fn_eeg,'file');
else MEG& EEG
    test = exist(fn_eeg,'file') & exist(fn_meg,'file');
end

if (OPTIONS.BEM.ForceXferComputation | ~test) | OPTIONS.BEM.Interpolative

    % Force (re)computation of transfer matrices, even if files exist in current study folder

    %    if OPTIONS.Verbose, bst_message_window('Computing the BEM Transfer Matrix (this may take a while)....'), end

    global nfv
    nfv = bem_xfer(R_eeg,R_meg1,O_meg1,Vertices,Faces,Param(1).Conductivity,Param(1).mode, ...
        Param(1).basis_opt,Param(1).test_opt,Param(1).ISA,fn_eeg,fn_meg,Param.Ntess_max,OPTIONS.Verbose,OPTIONS.BEM.checksurf);

    if ~isempty(find(flaggrad))
        if OPTIONS.Verbose, bst_message_window({'Gradiometers detected',...
                'Computing corresponding Gain Matrix. . .'}), end

        %fn_meg_2 = fullfile(pwd,[OPTIONS.rooot,'_megxfer_2.mat']);
        fn_meg_2 = sprintf('%s_megxfer2_%s_%s.mat',OPTIONS.rooot, OPTIONS.BEM.Basis,OPTIONS.BEM.Test);
        bem_xfer(R_eeg,R_meg2,O_meg2,Vertices,Faces,Param(1).Conductivity,Param(1).mode, ...
            Param(1).basis_opt,Param(1).test_opt,Param(1).ISA,fn_eeg,fn_meg_2,Param.Ntess_max,0,OPTIONS.BEM.checksurf); % Verbose = 0
        if OPTIONS.Verbose, bst_message_window('Gradiometer Channel Gain Matrix is Completed.'), end
    end

    if ~isempty(irefsens) % Do the same for reference channels
        %fn_meg_REF = fullfile(pwd,[OPTIONS.rooot,'_megxfer_REF.mat']);
        fn_meg_REF = sprintf('%s_megREFxfer_%s_%s.mat',OPTIONS.rooot, OPTIONS.BEM.Basis,OPTIONS.BEM.Test);
        bem_xfer(R_eeg,R_meg_REF,O_meg_REF,Vertices,Faces,Param(1).Conductivity,Param(1).mode, ...
            Param(1).basis_opt,Param(1).test_opt,Param(1).ISA,fn_eeg,fn_meg_REF,Param.Ntess_max,0,OPTIONS.BEM.checksurf);% Verbose = 0

        if ~isempty(find(flaggrad_REF))

            if OPTIONS.Verbose, bst_message_window({'MEG Reference Channels detected',...
                    'Computing corresponding Gain Matrix. . .'}), end

            %fn_meg_REF2 = fullfile(pwd,[OPTIONS.rooot,'_megxfer_REF2.mat']);
            fn_meg_REF2 = sprintf('%s_megREFxfer2_%s_%s.mat',OPTIONS.rooot, OPTIONS.BEM.Basis,OPTIONS.BEM.Test);

            bem_xfer(R_eeg,R_meg_REF2,O_meg_REF2,Vertices,Faces,Param(1).Conductivity,Param(1).mode, ...
                Param(1).basis_opt,Param(1).test_opt,Param(1).ISA,fn_eeg,fn_meg_REF2,Param.Ntess_max,0,OPTIONS.BEM.checksurf);% Verbose = 0
            if OPTIONS.Verbose, bst_message_window('MEG Reference Channel Gain Matrix is Completed.'), end
        end
    end
end

% Computation of TRANSFER MATRIX Completed ___________________________________________________


%%%% THIS PART SPECIFIES PARAMETERS USED TO GENERATE THE 3-D GRID %%%%%%%%%%


if OPTIONS.BEM.Interpolative%1%~exist(BEMGridFileName,'file')

    if OPTIONS.Verbose, bst_message_window('Computing BEM Interpolative Grid. . .'), end

    BEMGridFileName = [OPTIONS.rooot,'grid.mat'];

    % update tessellated envelope with surfaces possibly modified by
    % BEM_XFER (downsampling, alignment, embedding etc.)
    Vertices = {nfv(:).vertices};
    Faces = {nfv(:).faces};
    gridmaker(Vertices,Faces,BEMGridFileName,OPTIONS.Verbose);

    if OPTIONS.Verbose,
        bst_message_window('Computing BEM Interpolative Grid -> DONE'),

        % Visualization of surfaces + grid points
        for k = 1:length(Vertices)
            [hf,hs(k),hl] = view_surface('Head envelopes & a subset of BEM interpolative grid points',Faces{k},Vertices{k});
            view(90,0)
            delete(hl)
        end
        camlight
        rotate3d on

        set(hs(1),'FaceAlpha',.3,'edgealpha',.3,'edgecolor','none','facecolor','r')
        set(hs(2),'FaceAlpha',.2,'edgealpha',.2,'edgecolor','none','facecolor','g')
        set(hs(3),'FaceAlpha',.1,'edgealpha',.1,'edgecolor','none','facecolor','b')

        hold on

        load(BEMGridFileName)
        hgrid = scatter3(Rq_bemgrid(1:10:end,1),Rq_bemgrid(1:10:end,2),Rq_bemgrid(1:10:end,3),'.','filled');

    end

else

    global GBEM_grid

    % Grid points are the locations of the distributed sources
    if isempty(OPTIONS.Cortex) % CBB (SB, 07-May-2004)| Should work also for volumic grid
        errordlg('Please select a cortical grid for computation of BEM vector fields')
        return
    end

    try
        load(OPTIONS.Cortex.FileName); % Load tessellation supporting the source locations and orientations
    catch
        Users = get_user_directory;
        load(fullfile(Users.SUBJECTS,OPTIONS.Cortex.FileName)); % Load tessellation supporting the source locations and orientations
    end

    BEMGridFileName.Loc = Vertices{OPTIONS.Cortex.iGrid}; clear Vertices

    if OPTIONS.ApplyGridOrient % Take cortex normals into account
        ptch = patch('Vertices',BEMGridFileName.Loc','Faces',Faces{OPTIONS.Cortex.iGrid},'Visible','off');
        set(get(ptch,'Parent'),'Visible','off')
        clear Faces
        BEMGridFileName.Orient = get(ptch,'VertexNormals')';
        delete(ptch);
        [nrm,BEMGridFileName.Orient] = colnorm(BEMGridFileName.Orient);
        % Now because some orientations may be ill-set to [0 0 0] with the get(*,'VertexNormals' command)
        % set these orientation to arbitrary [1 1 1]:
        izero = find(nrm == 0);clear nrm
        if ~isempty(izero)
            BEMGridFileName.Orient(:,izero) = repmat([1 1 1]'/norm([1 1 1]),1,length(izero));
        end
        clear izero
    else
        BEMGridFileName.Orient = [];
    end

    %if OPTIONS.Verbose, bst_message_window(sprintf('Loading BEM interpolative grid points from %s', BEMGridFileName)), end

end



% This part computes the gain matrices defined on precomputed grid------------------------------------------------------------

% Assign file names where to store the gain matrices
if MEG & ~EEG
    bem_xfer_mfname = {fn_meg};
    %BEMGaingridFname = [OPTIONS.rooot,'_meggain_grid.mat'];
    BEMGaingridFname = sprintf('%s_MEGGainGrid_%s_%s.mat',OPTIONS.rooot, OPTIONS.BEM.Basis,OPTIONS.BEM.Test);
    test = exist(BEMGaingridFname,'file');
elseif EEG & ~ MEG
    bem_xfer_mfname = {fn_eeg};
    %BEMGaingridFname = [OPTIONS.rooot,'_eeggain_grid.mat'];
    BEMGaingridFname = sprintf('%s_EEGGainGrid_%s_%s.mat',OPTIONS.rooot, OPTIONS.BEM.Basis,OPTIONS.BEM.Test);
    test = exist(BEMGaingridFname,'file');
elseif EEG & MEG
    bem_xfer_mfname = {fn_meg,fn_eeg};
    %     BEMGaingridFname.MEG = [OPTIONS.rooot,'_meggain_grid.mat'];
    %     BEMGaingridFname.EEG = [OPTIONS.rooot,'_eeggain_grid.mat'];
    BEMGaingridFname.MEG = sprintf('%s_MEGGainGrid_%s_%s.mat',OPTIONS.rooot, OPTIONS.BEM.Basis,OPTIONS.BEM.Test);
    BEMGaingridFname.EEG = sprintf('%s_EEGGainGrid_%s_%s.mat',OPTIONS.rooot, OPTIONS.BEM.Basis,OPTIONS.BEM.Test);
    test = exist(BEMGaingridFname.MEG,'file') & exist(BEMGaingridFname.EEG,'file');
end

if 1%OPTIONS.BEM.ForceXferComputation | ~test% Recompute gain matrices when transfer matrices have been recomputed just before

    if OPTIONS.Verbose
        if OPTIONS.BEM.Interpolative
            bst_message_window('Computing the BEM gain matrix for interpolative grid. . .'),
        else
            bst_message_window('Computing the BEM gain matrix for source grid. . .'),
        end
    end

    t0 = clock;

    if OPTIONS.Verbose,
        if MEG & EEG
            bst_message_window('for MEG and EEG channels. . .')
        elseif EEG
            bst_message_window('for EEG channels. . .')
        elseif MEG
            bst_message_window('for MEG channels. . .')
        end
    end


    if length(bem_xfer_mfname) == 1 % xor(MEG,EEG)
        bem_gain(BEMGridFileName,bem_xfer_mfname{1},Param(1).ISA,BEMGaingridFname, OPTIONS.Verbose);
    else
        % MEG gaingrid matrix
        bem_gain(BEMGridFileName,bem_xfer_mfname{1},Param(1).ISA,BEMGaingridFname.MEG, OPTIONS.Verbose);
        % EEG gaingrid matrix
        bem_gain(BEMGridFileName,bem_xfer_mfname{2},Param(1).ISA,BEMGaingridFname.EEG, OPTIONS.Verbose);
    end

    if OPTIONS.Verbose, bst_message_window('-> DONE'), end

    if ~isempty(find(flaggrad)) % Compute forward model on the second set of magnetometers from the gradiometers array

        bem_xfer_mfname = sprintf('%s_megxfer2_%s_%s.mat',OPTIONS.rooot, OPTIONS.BEM.Basis,OPTIONS.BEM.Test);
        bem_gaingrid_mfname_2 = sprintf('%s_MEGGainGrid2_%s_%s.mat',OPTIONS.rooot, OPTIONS.BEM.Basis,OPTIONS.BEM.Test);

        if OPTIONS.Verbose, bst_message_window('MEG - completing gradiometers. . .'), end

        bem_gain(BEMGridFileName,bem_xfer_mfname,Param(1).ISA,bem_gaingrid_mfname_2,OPTIONS.Verbose);

        if OPTIONS.Verbose, bst_message_window('-> DONE'), end

        if MEG & EEG
            G1 = load(BEMGaingridFname.MEG);
        else
            G1 = load(BEMGaingridFname);
        end
        G2 = load(bem_gaingrid_mfname_2);
        % Apply respective weights within the gradiodmeters
        meg_chans = [OPTIONS.Channel(MEGndx(find(flaggrad))).Weight];
        w1 = meg_chans(1:2:end);
        w2 = meg_chans(2:2:end);

        G1.GBEM_grid(find(flaggrad),:) = w1'*ones(1,size(G1.GBEM_grid,2)).*G1.GBEM_grid(find(flaggrad),:)...
            +  w2' * ones(1,size(G1.GBEM_grid,2)).*G2.GBEM_grid(find(flaggrad),:);
        clear G2

        % is there special reference channel considerations?
        if ~isempty(irefsens)

            % Forward model on all reference sensors
            bem_xfer_mfname_REF = sprintf('%s_megREFxfer_%s_%s.mat',OPTIONS.rooot, OPTIONS.BEM.Basis,OPTIONS.BEM.Test);
            bem_gaingrid_mfname_REF = sprintf('%s_MEGGainGrid_REF_%s_%s.mat',OPTIONS.rooot, OPTIONS.BEM.Basis,OPTIONS.BEM.Test);

            if OPTIONS.Verbose, bst_message_window('MEG reference channels. . .'), end
            bem_gain(BEMGridFileName,bem_xfer_mfname_REF,Param(1).ISA,bem_gaingrid_mfname_REF,OPTIONS.Verbose);
            if OPTIONS.Verbose, bst_message_window('-> DONE'), end

            if ~isempty(find(flaggrad_REF)) % Gradiometers are present in reference channels

                bem_xfer_mfname_REF2 = sprintf('%s_megREFxfer2_%s_%s.mat',OPTIONS.rooot, OPTIONS.BEM.Basis,OPTIONS.BEM.Test);
                bem_gaingrid_mfname_REF2 = sprintf('%s_MEGGainGrid2_REF_%s_%s.mat',OPTIONS.rooot, OPTIONS.BEM.Basis,OPTIONS.BEM.Test);

                if OPTIONS.Verbose, bst_message_window('MEG reference channels / completing gradiometers. . .'), end

                bem_gain(BEMGridFileName,bem_xfer_mfname_REF2,Param(1).ISA,bem_gaingrid_mfname_REF2,OPTIONS.Verbose);

                if OPTIONS.Verbose, bst_message_window('-> DONE'), end

                GR = load(bem_gaingrid_mfname_REF);
                GR2 = load(bem_gaingrid_mfname_REF2);

                meg_chans = [OPTIONS.Channel(irefsens(find(flaggrad_REF))).Weight];
                w1 = meg_chans(1:2:end);
                w2 = meg_chans(2:2:end);

                GR.GBEM_grid(find(flaggrad_REF),:) = w1'*ones(1,size(GR.GBEM_grid,2)).*GR.GBEM_grid(find(flaggrad_REF),:)...
                    +  w2' * ones(1,size(GR.GBEM_grid,2)).*GR2.GBEM_grid(find(flaggrad_REF),:);

                %GR.GBEM_grid(find(flaggrad_REF),:) = GR.GBEM_grid(find(flaggrad_REF),:) - GR2.GBEM_grid(find(flaggrad_REF),:);
                clear GR2;
            end

            % Apply nth-order gradient correction on good channels only
            %Weight by the current nth-order correction coefficients
            %G1.GBEM_grid = G1.GBEM_grid - Channel(MEGndx(1)).Gcoef(find(ChannelFlag(irefsens)>0),:)*GR.GBEM_grid;

            G1.GBEM_grid = G1.GBEM_grid - OPTIONS.Channel(MEGndx(1)).Comment*GR.GBEM_grid;

        end

        GBEM_grid = 1e-7*G1.GBEM_grid;

        if MEG & EEG
            save(BEMGaingridFname.MEG,'GBEM_grid','-append');
        else
            save(BEMGaingridFname,'GBEM_grid','-append');
        end

    end

    % Detect EEG reference - if none is specified in the .Comment or Type fields,
    % and if none was passed through OPTIONS.EEGRef
    % -> Apply average-referencing of the potentials by default.

    if EEG
        % EEG Reference Channel
        EEGREFndx = good_channel(OPTIONS.Channel,[],'EEG REF');
        
        if MEG & EEG
            load(BEMGaingridFname.EEG,'GBEM_grid')
        else
            load(BEMGaingridFname,'GBEM_grid')
        end
        
        
        if isempty(EEGREFndx)% AVERAGE REF
            GBEM_grid = GBEM_grid - repmat(mean(GBEM_grid),size(GBEM_grid,1),1);
        else
           % GBEM_grid = GBEM_grid(setdiff(EEGndx,EEGREFndx)-EEGndx(1)+1,:) - repmat(GBEM_grid(EEGREFndx-EEGndx(1)+1,:),size(GBEM_grid(setdiff(EEGndx,EEGREFndx)-EEGndx(1)+1,:),1),1);
           GBEM_grid = GBEM_grid(2:end,:) - repmat(GBEM_grid(1,:),length(EEGndx),1); % SB : EEG REF is stored as first sensor in GBEM_grid; see line 2226
        end
        
        if MEG 
            save(BEMGaingridFname.EEG,'GBEM_grid','-append')
        else
            save(BEMGaingridFname,'GBEM_grid','-append')    
        end
        
    end
    
    
    if MEG & EEG
        meg = load(BEMGaingridFname.MEG,'GBEM_grid');
        eeg = load(BEMGaingridFname.EEG,'GBEM_grid');
        GBEM_grid = zeros(length(OPTIONS.Channel),size(GBEM_grid,2));
        GBEM_grid(MEGndx,:)= meg.GBEM_grid; clear meg
        GBEM_grid(EEGndx,:)= eeg.GBEM_grid; clear eeg
        
    elseif MEG
        
        meg = load(BEMGaingridFname,'GBEM_grid');
        GBEM_grid = zeros(length(OPTIONS.Channel),size(GBEM_grid,2));
        GBEM_grid(MEGndx,:)= meg.GBEM_grid; clear meg
        
    elseif EEG
        
        eeg = load(BEMGaingridFname,'GBEM_grid');
        GBEM_grid = zeros(length(OPTIONS.Channel),size(GBEM_grid,2));
        EEGndx = OPTIONS.EEGndx;
        GBEM_grid(EEGndx,:)= eeg.GBEM_grid; clear eeg
        
        %clear GBEM_grid
        %             % Now save the combined MEG/EEG gaingrid matrix in a single file 
        %             MEGEEG_BEMGaingridFname = strrep(BEMGaingridFname.EEG,'eeg','meg_eeg');
        %             Gmeg = load(BEMGaingridFname.MEG,'GBEM_grid');
        %             GBEM_grid = NaN * zeros(length(OPTIONS.Channel),size(Gmeg.GBEM_grid,2));
        %             GBEM_grid(MEGndx,:) = Gmeg.GBEM_grid; clear Gmeg 
        %             eeg = load(BEMGaingridFname.EEG);
        %             
        %             save_fieldnames(eeg,MEGEEG_BEMGaingridFname);
        %             
        %             GBEM_grid(setdiff(EEGndx,EEGREFndx),:) = eeg.GBEM_grid; clear eeg
        %             
        %             save(MEGEEG_BEMGaingridFname,'GBEM_grid','-append')
        %             
        %             BEMGaingridFname = MEGEEG_BEMGaingridFname;
    else
        save(BEMGaingridFname,'GBEM_grid','-append')
    end
    
    telap_meg_interp = etime(clock,t0);
    
    if OPTIONS.Verbose, bst_message_window(sprintf('Completed in %3.1f seconds', telap_meg_interp),...
            'Computing the Gain Matrix for the Interpolative Grid Points -> DONE'), end
else
    %     if MEG & EEG 
    %         if OPTIONS.Verbose, bst_message_window(sprintf('Loading 3D grid gain matrix from %s', BEMGaingridFname.EEG)),end 
    %     else
    %         if OPTIONS.Verbose, bst_message_window(sprintf('Loading 3D grid gain matrix from %s', BEMGaingridFname)),end 
    %     end
end



function g = gterm_constant(r,rq)
%gterm_constant 
% function g = gterm_constant(r,rq)

%<autobegin> -------- 20-Nov-2002 14:06:02 ------------------------------
% ---- Automatically Generated Comments Block using auto_comments -----------
%
% Alphabetical list of external functions (non-Matlab):
%   toolbox\rownorm.m
%<autoend> ---------- 20-Nov-2002 14:06:02 ------------------------------


if size(rq,1)  == 1 % Just one dipole
    r_rq= [r(:,1)-rq(1),r(:,2)-rq(2),r(:,3)-rq(3)];
    n = rownorm(r_rq).^3;
    g = r_rq./[n,n,n];
else
    g = zeros(size(r,1),3*size(rq,1));
    isrc = 1;
    for k = 1:size(rq,1)
       r_rq= [r(:,1)-rq(k,1),r(:,2)-rq(k,2),r(:,3)-rq(k,3)];
       n = rownorm(r_rq).^3;
       g(:,3*(isrc-1)+1: 3*isrc) = r_rq./[n,n,n];
       isrc = isrc + 1;
   end
   
end
