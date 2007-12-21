function varargout = spm_eeg_inv_model(action,varargin)

% spm_eeg_inv_model is a multi-purpose routine that deals with the generation
% of the head model for the solution of the forward problem.
%
% Called without arguments :
%
% >> model = spm_eeg_inv_model
% 
% the function calls GUI to select the operation(s) to be performed and
% to enter/select the required inputs.
% Here are the main options:
% - Generate scalp/brain vol & tessalate
% - Tessalate predefined binarized volume
% - Project electrodes coord (sphere or real) on scalp/brain surface
% - Define realistic sphere model
% - The 4 steps at once
%
% A 'model' structure is returned. Depending on the operation executed,
% the model will contain: surface meshes, electrodes information,
% realistic sphere model.
% 
% Most subfunctions can be called individually. Some even have their own
% little GUI because I used to use them on their own...
%
% The various structure formats (model, head, electrodes) are described 
% at the end of this help.
%
%__________________________________________________________________________
% FORMAT [Pvol,Pinv_sn,flags] = spm_eeg_inv_model('GenBin',Pvol,flag_bi);
%
% Generate the binarised image of the brain (and scalp if the image allows).
% 
% Input :
%   Pvol    : filename of the image used to build the head model.
%   flag_bi : useful flags
%      * img_type: image type 1=T1 MRI, 2 = PET, 3 = EPI
%      * img_norm: flag indicating if the image is normalised (1) or not (0)
%      * ne,ng   : # iteration for erosion and growing
%      * thr_im  : threhold applied to binarise image
%   NOTE: these flags are not always necessary anymore...
% Output :
%   Pvol    : file name of the bin images gnerated, scalp and iskull
%   Pinv_sn : file name of mat file containing inverse of the normalisation
%             transform.
%   flags   : all the flags with added file names for original image, brain
%             mask and spatial transformation
%__________________________________________________________________________
%
% FORMAT [head,flags] = spm_eeg_inv_model('GenMesh',Pvol,flags);
%
% Generate the tesselated surface of the 1 to 3 head volumes
% (brain-skull-scalp) from the binarized volumes.
% Input :
%   Pvol    : filenames of 3 head volumes 'brain', 'skull' & 'scalp'.
%   flags    : various flags
%      * n       : (1x3) provides the number of vertices on each surface
%                   Npt = n(i)^2*5/4+2
%      * br_only : only use the binarised brain volume (1) or scalp (0, default)
%      * q_elastm: Correct mesh using an elastic model (1, default), or not (0).
%      * q_meased: Measure edges of mesh (1), or not (0, default)
%      * q_4thpt : Determine 4th point central of each traingle (1), or not (0, default)
% Output :
%   head    : head structures with 1 to 3 tesselated surfaces
%   flags   : all the flags used
%__________________________________________________________________________
%
% FORMAT [electrodes,flags] = spm_eeg_inv_model('Elec2Scalp',surf,el_loc,el_name,flags);
%
% Projecting the electrodes (in spherical or realistic coordinates) 
% on the realistic head model.
% 
% Input :
% - surf    : surface on which the electrodes are projected, 
%             if it's the brain surface, then the scalp-brain width is added
% - el_loc  : provided electrode coordinates, on a standard sphere or
%             real coordinates as measured estimated on the patient's head.
% - el_name : electrode names.
% - flags   : various option flags and extra bits of information
%   * q_RealLoc : 'el_loc' contains coordinates on a standard sphere (0),
%                 or as measured on the patient's head (1)
%   * q_RealNI  : nasion/inion coordinates are from the template (0),
%                 or as estimated on the patient's image
%   * br_only   : only the brain surface is available (1),
%                 or the scalp surface can be used (0)
%   * nasion    : coordiantes (in mm), on the template (defaults value filled
%   * inion     : in if left empty) or as estimated on the patient's image
%   * scbr_w    : scalp-brain surface width, in mm. (Default 20)
%   * Mtempl    : Affine transform mapping from the patient image space
%                 into the template space
% Output :
%   electrodes  : electrode structure
%   flags       : all the flags used
%
% When dealing with standard electrode locations on a sphere, one uses 
% the coordinates of the nasion and inion to estimate the pitch angle, then
% the electrode locations are mapped into the patient space
% When the real (approximate) electrode locations is provided, these are 
% simply adjusted on the scalp surface
% If the model is 'brain only' (from a PET or EPI scan), the approximate 
% brain-scalp width is required.
%
% The electrodes substructure is created at the end.
%__________________________________________________________________________
%
% FORMAT [pt_mm] = spm_eeg_inv_model('vx2mm',pt_vx,M);
%
% Transforms point(s) 'pt_vx' voxel coordinates into mm coordinates 'pt_mm'
% according to the transformation matrix 'M'.
% 'M' is the 4x4 affine transformation matrix : from vx to mm
%__________________________________________________________________________
%
% FORMAT [pt_vx] = spm_eeg_inv_model('mm2vx',pt_mm,M);
%
% Transforms point(s) 'pt_mm' mm coordinates into voxel coordinates 'pt_vx'
% according to the transformation matrix 'M'.
% 'M' is the 4x4 affine transformation matrix : from vx to mm
%__________________________________________________________________________
%
% FORMAT [Centre,c_mm] = spm_eeg_inv_model('CtrBin',P);
%
% Determine the "centre of mass" (vx and mm) of a bin image.
%__________________________________________________________________________
%
% FORMAT [vn,vn_i,S,S_i] = spm_eeg_inv_model('NormTri',l_tr,ts)
%
% Calcultates the normal and surface of a triangle on a tesselated surface.
% If a list of triangle is provided, then the "mean" normal and the normals are calculated.
% Input :
%     l_tr    = index of triangle
%     ts   = tesselated surface
% Output:
%     vn   = (mean) normalised normal vector
%     vn_i = individual normal vectors
%     S    = (mean) surface
%     S_i  = individual surfaces
%__________________________________________________________________________
%
% FORMAT ts = spm_eeg_inv_model('TesBin',n,Centre,P,info);
%
% Generate a mesh covering a binarized volume
%   1. Generate a simple spherical mesh
%   2. The spherical mesh is projected radially on the bin volume
% Afterwards, "elastic" mesh correction can thus be useful to correct 
% some overlong edges.
%
% Input : 
%   n:      number of vertices on each surface Npt = n^2*5/4+2
%   Centre: centre of bin volume for the radial projection
%   P:      filename of bin volume
%   info:   information string
%__________________________________________________________________________
%
% FORMAT ts = spm_eeg_inv_model('ElastM',ts);
%
% Modify the the mesh in order to reduce overlong edges.
% The procedure uses an elastic model :
% At each vertex, the neighbouring triangles and vertices connected directly
% are considered.
% Each edge is considered elastic and can be lengthened or shortened,
% depending on their legnth.
% Algorithm: G.Taubin, A signal processing approach to fair surface design, 1995
% This is a non-shrinking smoothing algo.
%
% Input : 
%   ts:     tesselated surface
% Output :
%   ts:     tesselated surface with corrected mesh
%__________________________________________________________________________
%
% FORMAT ts = spm_eeg_inv_model('MeasEd',ts);
%
% This measures the edges of the mesh.
%   Input : ts
%   Output : ts with measured edges
%__________________________________________________________________________
%
% FORMAT ts = spm_eeg_inv_model('Tr4thPt',ts);
%
% Finds the 4th' point of each triangle:
% a flat triangle is only an approximation of the real underlying volume surface, 
% a 4th point taken at the centre of the triangle is found by approximating
% the neighbouring triangles in a sphere, to provide an idea of the local curvature.
% This is used in the BEM solution of the forward problem, to sort out 
% the problem of the auto-solid angle.
%
% Input : 
%   ts:     tesselated surface
% Output :
%   ts:     tesselated surface with 4th point at each triangle
%__________________________________________________________________________
%
% FORMAT tsph = spm_eeg_inv_model('TesSph',n,r);
%
% Generate a structure 'tsph' containing a tesselated sphere.
% Input : 
%   n : number of 'latitude' divisions on the sphere. It MUST be even!
%   r : radius of the sphere
% Output : 
%   tsph .vert : vertices coordinates (3 x Nvert)
%		 .tri  : triangle patches (3 x Ntri)
%		 .info : info string
%__________________________________________________________________________
%
% FORMAT [Pout/val] = spm_eeg_inv_model('ErodeGrow',P/val,ne,ng,thr_im)
% 
% It erodes then grows an image after thresholding it.
% Inputs : P/val,ne,ng,thr_im
%     P : file name of image to erode/grow
%     val : full volume of image (as loaded by spm_read_vols)
%     ne : nr of steps for erosion, default 3
%     ng : nr of steps for growing, default 6
%     thr_im : threshold value to apply, default .8
%               (compared to the maximal value of image)
% Output : Pout/val
%     Pout : name of the file generated
%     val : full value of eroded-grown image
%     
% Notes:
% - if 2 file names are specified in P[input], the 2nd one is used as the name
%   of the generated file otherwise the new filename is creted from P as
%           [P,'_e',num2str(ne),'g',num2str(ng),'.img']
% - If a file name is passed the output is a filename.
%   If a matrix of values is passed, the output is a matrix of values.
%
% IMPORTATNT !!!
% - It uses brutal force as it loads the whole image into memory !!!
%   things are eased up by usint unsigned integers.
% - voxels are supposed to be isotropic in size...
%
%_______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Christophe Phillips,
% $Id: spm_eeg_inv_model.m 1039 2007-12-21 20:20:38Z karl $

% Format of 'model' structure :
% #############################
%
% model.
%	head		: structure of head surfaces (1x3)
%	electrodes	: structure for electrodes
%   Mtempl      : affine transformation from the patient into the template space
%	sigma		: conductivity values (1x3)
%	(IFS		: cells containing Intermediate Forward Solution (1x3),
%			      these do NOT depend on the dipoles grid.)
%                 Only used for the BEM solution !
%
% 'head' sub-structure format
% ===========================
% .head(i)		: i=1,2,3 for brain, skull & scalp surfaces
%	XYZmm		: coordinates of vertices (3xNv)
%	tri		: indices of vertices in triangles (3xNt)
%	M		: voxel-to-mm space transformation matrix
%	nr		: [nr_of_vertices (Nv) & nr_of_triangles (Nt)]
%	pt4.XYZvx	: coord of projected triangle centre (3xNt)
%	edges.		: triangles edges information structure
%		ed	: edges length, vertices indices and triangles indices
%		nr	: [nr_edges mean_length std_length]
%	info		: informations
%
%
%	Format of .edges :
%	------------------
%	.edges.ed :					
% 		ed(:,1) = length,
%		ed(:,[2 3]) = indices of 2 vertices, e.g ind. of a and c	
%		ed(:,[4 5]) = indices of 2 triangles, e.g ind. of I and II
%
%		  length	    a x------x b
%		 _______	     /\  I  /
%		'       `	    /  o   /
%		x---o---x	   / II \ /
%				    d x------x c
%
%	.edges.nr :
%		nr(1) : number of edges,
%		nr([2 3]) : mean and std of edges length.
%
% 'electrodes' sub-structure format
% =================================
% .electrodes.		: structure for electrodes
%	vert	: index of the closest vertex (Nel x 1)
%	tri		: index of the triangle beneath each electrode
%	nr		: nr of electrodes
%	XYZmm	: coord. of electrodes on the scalp (in voxel)
%	M		: voxel-to-mm space transformation matrix
%	Cnames	: electrodes names
%	info	:
%   nasion  :
%   inion   :

spm('FigName','Realistic head model');

if nargin == 0, action = 'Init'; end;

switch lower(action),

%________________________________________________________________________
case 'init'
%------------------------------------------------------------------------
% Use a GUI.
%------------------------------------------------------------------------
	pos = 1 ;

    % Select the operation through the GUI or from passed variable
    if nargin==2 && isnumeric(varargin{1})
        gener = round(varargin{1});
        if gener<1 || gener>5
            gener = 5;
            warning('Wrong options for model, I assume "generate all at once"');
        end
    else
		text_opt = [...
            'Generate scalp/brain vol & tessalate|',...
            'Tessalate predefined binarized volume|',...
			'Project electrodes coord (spher or real) on scalp/brain surface|',...
			'Define realistic sphere model|',...
            'The 4 steps at once'];
		gener = spm_input('Generate what ?',pos,'m',text_opt);
    end

    % Load required inputs
    %======================
    if gener==1 || gener==5
    % Generate scalp/brain vol & tessalate
        Pvol = spm_select(1,'image','Image to build model');
%         flag_bi.img_type = spm_input('Image type :','+1','T1 MRI|PET|EPI',[1,2,3],1);
        flag_bi.img_type = 1; % Assume it's always a structural MR image.
%         flag_bi.img_norm = spm_input('Image normalised already :','+1','y/n',[1,0],2);
        flag_bi.img_norm = 0; % Assume it is not normalised yet.
        if flag_bi.img_type==1
            flag_te.br_only=0;
            flag_te.Nvol = 2;
        else
            flag_te.br_only=1;
            flag_te.Nvol = 1;
        end
    end
        
    if gener==2
    % Tessalate predefined binarized volume only.
		Pvol = spm_select(Inf,'*_o*.img','1-2-3 bin images to tessalate (br-sk-sc)');
        flag_te.Nvol = size(Pvol,1);
        if flag_te.Nvol==1
            flag_te.br_only = spm_input('Is this the brain bin volume ?','+1','y/n',[1,0],0);
        else
            flag_te.br_only=0;
        end
        if ~isempty(PMtempl)
            load(PMtempl)
        else
            warning('It is assumed that the bin volumes are in the template space !')
            Mtempl = eye(4);
        end
    end
    if gener==1 || gener==2 || gener==5
    % Generate scalp/brain vol & tessalate binarized volume
        if flag_te.Nvol>1 || flag_te.br_only
%     		Npt = spm_input('Approx. nr of vertices on brain surf','+1','e',4000);
    		Npt = 4000; % 4000 seems a "good" number for the inner skull surface
            flag_te.n = zeros(1,3); flag_te.n(1) = 2*round((sqrt((Npt-2)*4/5))/2);
        end
        if ~flag_te.br_only 
%     		Npt = spm_input('Approx. nr of vertices on sk/sc surf','+1','e',2000);
    		Npt = 4000; % Use also about 4000 vertices for the scalp surface
            if flag_te.Nvol>1 
                flag_te.n(2:3) = 2*round((sqrt((Npt-2)*4/5))/2);
            else
                flag_te.n = 2*round((sqrt((Npt-2)*4/5))/2);
            end
        end
%         flag_te.q_elast_m  = spm_input('Correct the projected mesh ?','+1','y/n',[1,0],1) ;
        flag_te.q_elast_m  = 1 ; % of course do the "surface smoothing".
        [pth,fn,ext,nr] = spm_fileparts(Pvol);
        Pmod = fullfile(pth,['model_head_',fn,'.mat']);
	end

	if gener==3 || gener==5
    % Project electrodes coord (spher or real) on scalp/brain surface
        if gener==3
            Pmod = spm_select(1,'^model.*\.mat$','Head model file');
            load(Pmod)
        end
        if exist('flag_te') && isfield(flag_te,'br_only')
            flags_el.br_only = flag_te.br_only;
        elseif exist('model') && isfield(model.param,'br_only');
            flags_el.br_only = model.param.br_only;
        else
            if length(model.head)>1
                flags_el.br_only = 0;
            else
                flags_el.br_only = spm_input('Only one surface modeled :','+1',...
                                'scalp|brain',[0,1],1);
                if flags_el.br_only
                    flags_el.scbr_w = spm_input('Scalp-brain width :','+1','e',20);
                end
            end
            model.br_only = flags_el.br_only;
        end
        if exist('Mtempl')
            flags_el.Mtempl = Mtempl;
        elseif exist('model') && isfield(model.param,'Mtempl');
            flags_el.Mtempl = model.param.Mtempl;
        else
            warning('It is assumed that your original image was in template space.')
        end

        q_StdElec = spm_input('Electrode set :','+1','standard|user defined',[1,0],1);
        if q_StdElec
            % Select standard electrode set, as defined in 'spm_eeg_inv_electrset'
            [set_Nel,set_name] = spm_eeg_inv_electrset;
			el_set = spm_input('Which set of electrodes', ...
					'+1','m',set_name);
     		[el_sphc,el_name] = spm_eeg_inv_electrset(el_set) ;
            flags_el.q_RealLoc = 0;
        else
            % Enter user defined electrode set.
            el_sphc = spm_input('All electrode coordinates','+1','e')
            el_name = spm_input('Electrod names (cell array)','+1','e')
            flags_el.q_RealLoc = spm_input('Electrode set :','+1', ...
                            'real location|on a sphere',[1,0],1);
            flags_el.q_RealNI  = spm_input('Do you have the nasion/inion location',...
                            '+1','y/n',[1,0],2);
            if flags_el.q_RealNI
                flags_el.nasion = spm_input('Nasion coordinates (in mm)','+1','e',[],3);
                flags_el.inion  = spm_input('Inion coordinates (in mm)','+1','e',[],3);
            end
        end
	end

    if gener==4
    % Define realistic sphere
        Pmod = spm_select(1,'^model.*\.mat$','Head model file');
    end


    % Create things as selected
    %==========================
    if gener == 1 || gener == 5
        % Segment image, and create brain volume
        [Pbin,Pinv_sn,flags_bi] = spm_eeg_inv_model('GenBin',Pvol,flag_bi);
        if gener == 1 
            varargout{1} = Pbin;
            varargout{2} = Pinv_sn;
            varargout{3} = flags_bi;
        end
    end
    
    if gener==1 || gener==2 || gener==5
    % Tessalate predefined binarized volume
        [head,cc,flags_te] = spm_eeg_inv_model('GenMesh',Pbin,flag_te); % n,br_only,q_elast_m
        model.head = head;
        if exist('flags_bi'), model.flags.fl_bin = flags_bi; end
        model.flags.fl_tess = flags_te;
        model.param = struct('br_only',0,'Nvol',1,'Mtempl',[], ...
                             'Pinv_sn','','sigma',[.33 .004 .33]);
        if exist('Mtempl'), model.param.Mtempl = Mtempl; end
        if isfield(flag_te,'br_only'), model.param.br_only = flag_te.br_only; end
        if isfield(flag_te,'Nvol'), model.param.Nvol = flag_te.Nvol; end
        model.fname = Pmod;
        save(Pmod,'model')
        if gener==2 || gener==5
            varargout{1} = model ;
        elseif gener==1
            varargout{4} = model ;
        end
    end
        
    if gener==3 || gener==5
    % Project electrodes coord (spher or real) on scalp/brain surface
		fname_el = ['el_sphc',num2str(size(el_sphc,2))];
		save(fname_el,'el_sphc','el_name');
        [electrodes,flags_el] = spm_eeg_inv_model('Elec2Scalp',model.head(end), ...
                                                el_sphc,el_name,flags_el);
        model.electrodes = electrodes ;
        model.flags.fl_elec = flags_el ;
		Pmod_s = [spm_str_manip(Pmod,'s'),'_e', ...
				num2str(model.electrodes.nr),'.mat'] ;
		save(Pmod_s,'model')
        varargout{1} = model ;    
    end
    
    if gener==4 || gener==5
    % Define realistic sphere model
        if gener==4, load(Pmod); end
        if isfield(model,'br_only')
            flags_Rs.br_only = model.flags.fl_tess.br_only;
        else
            flags_Rs.br_only = 0;
        end
        [spheres,a,b,c,flags_Rs] = spm_eeg_inv_Rsph(model,[],flags_Rs);
        model.spheres = spheres;
        model.flags.fl_RealS = flags_Rs;
		Pmod_Rs = [spm_str_manip(Pmod,'s'),'_Rs.mat'] ;
		save(Pmod_Rs,'model')
        varargout{1} = model ;  
    end


%________________________________________________________________________
case 'genbin'
%------------------------------------------------------------------------
% FORMAT [Pbin,Pinv_sn,flags] = spm_eeg_inv_model('GenBin',Pvol,flag_bi);
%
% Generate the binarised image of the brain (and scalp if the image allows).
% 
% Input :
%   Pvol    : filename of the image used to build the head model.
%   flag_bi : useful flags
%      * img_type: image type 1=T1 MRI, 2 = PET, 3 = EPI
%      * img_norm: flag indicating if the image is normalised (1) or not (0)
%      * ne,ng   : # iteration for erosion and growing
%      * thr_im  : threhold applied to binarise image
%   NOTE: these flags are not always necessary anymore...
% Output :
%   Pbin    : file name of the bin images generated, scalp and iskull
%   Pinv_sn : file name of mat file containing inverse of the normalisation
%             transform.
%   flags   : all the flags with added file names for original image, brain
%             mask and spatial transformation
%------------------------------------------------------------------------
    fprintf('\n\n'), fprintf('%c','='*ones(1,80)), fprintf('\n')
    fprintf(['\tGenerating binarized volumes from image: Please wait.\n']);
    Pvol = varargin{1} ;
    def_flags = struct('img_type',1,'img_norm',0,'ne',1,'ng',2,'thr_im',[.5 .1]);
    if nargin<3
        flags = def_flags;
    else
        flags = varargin{2};
        fnms  = fieldnames(def_flags);
		for i=1:length(fnms),
			if ~isfield(flags,fnms{i}),
                flags = setfield(flags,fnms{i},getfield(def_flags,fnms{i}));
            end
		end
    end
    
    load('defaults_eeg_mesh.mat');

% 1.Segment image, and estimate the mapping between template and image space.
    jobs{1}.spatial{1}.preproc.data = Pvol;
    jobs{1}.spatial{1}.preproc.output.biascor = 1;
    jobs{1}.spatial{1}.preproc.output.GM = [0 0 1];
    jobs{1}.spatial{1}.preproc.output.WM = [0 0 1];
    jobs{1}.spatial{1}.preproc.output.CSF = [0 0 1];
    jobs{1}.spatial{1}.preproc.output.cleanup = 0;

    res = spm_preproc(Pvol);
    [sn,isn] = spm_prep2sn(res);
    [pth,nam,ext] = spm_fileparts(Pvol);
    Pinv_isn = fullfile(pth,[nam '_vbm_inv_sn.mat']);
    Pinv_sn = fullfile(pth,[nam '_vbm_sn.mat']);
    savefields(Pinv_isn,isn);
    savefields(Pinv_sn,sn);
    spm_preproc_write(sn,jobs{1}.spatial{1}.preproc.output);
        % write the segmented images in native space + bias corrected MRI

% 2.Adds up GM/WM/CSF to produce inner skull surface
    Pin   = strvcat(fullfile(pth,['c1',nam,ext]), ...
                    fullfile(pth,['c2',nam,ext]), ...
                    fullfile(pth,['c3',nam,ext]));
    Pisk  = fullfile(pth,[nam,'_iskull',ext]);
    fl_ic = {[],[],'uint8',[]};
    Pisk  = spm_imcalc_ui(Pin,Pisk,'i1+i2+i3',fl_ic);
    
% 3.Use a little bit of erosion/growing to refine the model
%   and write the *_iskull img on disk
    Parg = strvcat(Pisk,Pisk);
    ne   = flags.ne(1); ng = flags.ng(1); thr_im = flags.thr_im(1);
	[Pout] = spm_eeg_inv_model('ErodeGrow',Parg,ne,ng,thr_im);

    if flags.img_type==1
% 4.Generate the outer-scalp volume, if possible
        Pvolc   = fullfile(pth,['m',nam,ext]);
		Vsc     = spm_vol(Pvolc);
		Vsc.dat = spm_loaduint8(Vsc);
%         [mnv,mxv] = spm_minmax(Vsc.dat)
		ne = flags.ne(end); ng = flags.ng(end); thr_im = flags.thr_im(end);
		Vsc.dat = spm_eeg_inv_model('ErodeGrow',Vsc.dat,ne,ng,thr_im);
		
        % The bottom part of the image needs to be masked out, in MNI space.
% 		p1 = [0 85 -50]'; % point below on the nose
% 		p2 = [0 -55 -70]'; % point below the bottom of cerebellum (mm)
		p1 = [0 110 -105]'; % point below on the nose
		p2 = [0 -55 -110]'; % point below the bottom of cerebellum (mm)
		c = [p2(1:2)' ((p1(2)-p2(2))^2+p1(3)^2-p2(3)^2)/(2*(-p2(3)+p1(3)))];
		R2 = (p2(3)-c(3))^2;
		% center and radius of sphere to chop off bottom of the head
		
		X = (1:Vsc.dim(1))'*ones(1,Vsc.dim(2)); X =X(:)';
		Y = ones(Vsc.dim(1),1)*(1:Vsc.dim(2)); Y = Y(:)';
		Z = zeros(Vsc.dim(1),Vsc.dim(2)); Z = Z(:)'; 
		Unit = ones(size(Z));
		% X,Y,X coordinates in vox in original image
        isn_tmp = isn; isn_tmp.Tr = []; % using the simple affine transform should be enough

        for pp = 1:Vsc.dim(3)
            XYZ = Vsc.mat*[X ; Y ; Z+pp ; Unit]; % coord in mm of voxels of scalp img
            XYZ = spm_get_orig_coord(XYZ(1:3,:)', isn_tmp)';           
            lz = find(((XYZ(1,:)-c(1)).^2+(XYZ(2,:)-c(2)).^2+(XYZ(3,:)-c(3)).^2-R2)>0);
            val_pp = Vsc.dat(:,:,pp);
            val_pp(lz) = 0;
            Vsc.dat(:,:,pp) = val_pp;
		end
		
        % Write the file
        Vsc.fname = fullfile(pth,[nam,'_scalp',ext]);
        Vsc.dt = [2 0];
        Vsc.pinfo = [1/255 0 0]';
		Vsc = spm_create_vol(Vsc);
		spm_progress_bar('Init',Vsc.dim(3),'Writing Outer-scalp','planes completed');
		for pp=1:Vsc.dim(3),
			Vsc = spm_write_plane(Vsc,double(Vsc.dat(:,:,pp)),pp);
                % No need to divide by 255 as output from Erode grow is [0 1].
			spm_progress_bar('Set',pp);
		end;
		spm_progress_bar('Clear');
        varargout{1} = strvcat(Pisk,Vsc.fname);
    end
    Pinv_sn = strvcat(Pinv_isn,Pinv_sn);
    varargout{2} = Pinv_sn;
    flags.Pinv_sn = Pinv_sn;
    flags.Pimg = Pvol;
    flags.Pisk = Pisk
    varargout{3} = flags;
    
%________________________________________________________________________
case 'vx2mm'
%------------------------------------------------------------------------
% FORMAT [pt_mm] = spm_eeg_inv_model('vx2mm',pt_vx,M);
% Transforms point(s) 'pt_vx' voxel coordinates into mm coordinates 'pt_mm'
% according to the transformation matrix 'M'.
% 'M' is the 4x4 affine transformation matrix : from vx to mm
%------------------------------------------------------------------------
	pt_vx = varargin{1} ;
	M = varargin{2} ;
	
	if size(pt_vx,1)==3
        Npt = size(pt_vx,2);
	elseif size(pt_vx,2)==3;
        pt_vx = pt_vx';
        Npt = size(pt_vx,2);
	else
        error('Wrong vectors format !')
	end
	
	pt_mm = M*[pt_vx;ones(1,Npt)];
	varargout{1} = pt_mm(1:3,:);
%________________________________________________________________________
case 'mm2vx'
%------------------------------------------------------------------------
% FORMAT [pt_vx] = spm_eeg_inv_model('mm2vx',pt_mm,M);
% Transforms point(s) 'pt_mm' mm coordinates into voxel coordinates 'pt_vx'
% according to the transformation matrix 'M'.
% 'M' is the 4x4 affine transformation matrix : from vx to mm
%------------------------------------------------------------------------
	pt_mm = varargin{1} ;
	M = varargin{2} ;
	
	if size(pt_mm,1)==3
        Npt = size(pt_mm,2);
	elseif size(pt_mm,2)==3;
        pt_mm = pt_mm';
        Npt = size(pt_mm,2);
	else
        error('Wrong vectors format !')
	end
	
	pt_vx = M\[pt_mm;ones(1,Npt)];
	varargout{1} = pt_vx(1:3,:);
%________________________________________________________________________
case 'ctrbin'
%------------------------------------------------------------------------
% FORMAT [Centre,c_mm]=spm_eeg_inv_model('CtrBin',P);
% Determine the "centre of mass" (vx and mm) of a bin image.
%------------------------------------------------------------------------
	if nargin<2
		P = spm_select(1,'o*.img','Bin image to use');
    else
        P = varargin{1};
	end
	Vp = spm_vol(P) ;

	Vp.dat = spm_loaduint8(Vp);
	X = (1:Vp.dim(1))'*ones(1,Vp.dim(2)); X =X(:)';
	Y = ones(Vp.dim(1),1)*(1:Vp.dim(2)); Y = Y(:)';
	Z = zeros(Vp.dim(1),Vp.dim(2)); Z = Z(:)'; 
	Unit = ones(size(Z));
	
	c_mm = zeros(1,3); SI = 0;
	for pp=1:Vp.dim(3)
        I_pl = double(Vp.dat(:,:,pp)); I_pl = I_pl(:)*ones(1,3);
        XYZ = Vp.mat*[X ; Y ; Z+pp ; Unit];
        c_mm = c_mm + sum(XYZ(1:3,:)'.*I_pl);
        SI = SI+sum(I_pl(:,1));
	end
	c_mm = c_mm/SI;
    Centre = spm_eeg_inv_model('mm2vx',c_mm,Vp.mat)';
    varargout{1} = Centre;
    varargout{2} = c_mm;
    
%________________________________________________________________________
case 'elec2scalp'
%------------------------------------------------------------------------
% FORMAT [electrodes,flags] = spm_eeg_inv_model('Elec2Scalp',surf,el_loc,el_name,flags);
%
% Projecting the electrodes (in spherical or realistic coordinates) 
% on the realistic head model.
% 
% IN :
% - surf    : surface on which the electrodes are projected, 
%             if its the brain surface, then the scalp-brain width is added
% - el_loc  : provided electrode coordinates, on a standard sphere or
%             real coordinates as measured estimated on the patient's head.
% - el_name : electrode names.
% - flags   : various option flags and extra bits of information
%   * q_RealLoc : 'el_loc' contains coordinates on a standard sphere (0),
%                 or as measured on the patient's head (1)
%   * q_RealNI  : nasion/inion coordinates are from the template (0),
%                 or as estimated on the patient's image
%   * br_only   : only the brain surface is available (1),
%                 or the scalp surface can be used (0)
%   * nasion    : coordiantes (in mm), on the template (defaults value filled
%   * inion     : in if left empty) or as estimated on the patient's image
%   * scbr_w    : scalp-brain surface width, in mm. (Default 20)
%   * Mtempl    : Affine transform mapping from the patient image space
%                 into the template space
%
% When dealing with standard electrode locations on a sphere, one uses 
% the coordinates of the nasion and inion to estimate the pitch angle, then
% the electrode locations are mapped into the patient space
% When the real (approximate) electrode locations is provided, these are 
% simply adjusted on the scalp surface
% If the model is 'brain only' (from a PET or EPI scan), the approximate 
% brain-scalp width is required.
%
% The electrodes substructure is created at the end.
%------------------------------------------------------------------------
	fprintf('\n'), fprintf('%c','='*ones(1,80)), fprintf('\n')
    fprintf(['Placing electrodes on the scalp surface.\n']);
    fprintf('%c','='*ones(1,80)), fprintf('\n')
    surf   = varargin{1};
	el_loc = varargin{2};
	Nel    = size(el_loc,2);
    if nargin<4
        [set_Nel,set_name] = spm_eeg_inv_electrset;
        try 
            [kk,el_name] = spm_eeg_inv_electrset(find(set_Nel==Nel));
        catch
            el_name = [];
        end
    else
        el_name = varargin{3};
    end

    def_flags = struct('q_RealLoc',1,'br_only',0,'nasion',[0 84 -28],'inion',[0 -118 -26], ...
                        'q_RealNI',0,'scbr_w',[16],'Mtempl',eye(4));
    if nargin<5
        flags = def_flags;
    else
        flags = varargin{4};
        fnms  = fieldnames(def_flags);
		for i=1:length(fnms),
			if ~isfield(flags,fnms{i}),
                flags = setfield(flags,fnms{i},getfield(def_flags,fnms{i}));
            end
		end
    end
    
    if ~flags.q_RealLoc
        % Standard electrode set location is used here
        % and they're provided on a simple sphere
        % => need to use nasion and inion to orient them properly.

        % Centre translation, scaling factor and pitch angle
	    % 	-> afine transformation, correcting for the pitch ONLY!
        [na_mm] = flags.nasion;
        [in_mm] = flags.inion;
    	cNI = (na_mm+in_mm)/2 ;
    	s = norm(cNI-na_mm) ;
    	theta = acos( (na_mm(2)-cNI(2))/s ) ; % pitch angle to apply
        if flags.q_RealNI
        	Tr = spm_matrix([cNI theta 0 0 s s s]) ; % apply pitch correction
        else
        	Tr = inv(flags.Mtempl) * spm_matrix([cNI theta 0 0 s s s]) ;
            % same but takes into account the mapping with template space.
        end            

		% Electrodes: "realistic" sphere coordinates...
		el_rsc = Tr * [el_loc ; ones(1,Nel)] ;
		el_rsc = el_rsc(1:3,:);
    else
        el_rsc = el_loc;
        cNI    = surf.Centre;
    end
    % figure, plot3(el_rsc(1,:)',el_rsc(2,:)',el_rsc(3,:)','*'), axis equal, xlabel('axe x'),ylabel('axe y')
    
	electrodes.vert = zeros(Nel,1);
	electrodes.tri  = zeros(Nel,1);
	pos_el_mm       = zeros(3,Nel);

	% projection of the el_rsc on the scalp surface.
	% The "closest" point is detemined by the angle between the electrode_i
	% and the scalp vertices from the centre
	% Then the supporting triangle is found, such that the electrode_i falls
	% within the traingle and the "exact" coord of electrode_i is obtained.
	XYZmm = surf.XYZmm ;
	v_sv = XYZmm-(cNI'*ones(1,size(XYZmm,2))) ;
	v_sv = v_sv./( ones(3,1)*sqrt(sum(v_sv.^2)) ) ; % direction of surf vertices from cNI
	for ii=1:Nel
		v_el = (el_rsc(:,ii)-cNI')/norm(el_rsc(:,ii)-cNI') ;
			% orient vect. of electrode_i
			% orient vect. of scalp vertices
		alpha = acos(v_el'*v_sv) ;
		[m_a,ind_v] = min(abs(alpha)) ; % Find the vertex with closest orientation
		list_t = find( (surf.tri(1,:)==ind_v) | ...
				(surf.tri(2,:)==ind_v) | ...
				(surf.tri(3,:)==ind_v) ) ; % Triangles linked to that vertex.
		Nlist = length(list_t) ;
		is_in = zeros(Nlist,1) ;
		proj_el = zeros(3,Nlist) ;
		for j=1:Nlist
			v1 = XYZmm(:,surf.tri(1,list_t(j))) ;
			v2 = XYZmm(:,surf.tri(2,list_t(j))) ;
			v3 = XYZmm(:,surf.tri(3,list_t(j))) ;
			v13 = v3-v1 ; v12 = v2-v1 ;
			A = [v13 v12 v_el] ;
			b = el_rsc(:,ii) - v1;
			x = A\b ;
			proj_el(:,j) = el_rsc(:,ii) - x(3)*v_el ;
            % In fact, el_rsc(:,ii) = v1 + v13*x(1) + v12*x(2) + v_el*x(3)
            % and proj_el(:,j) = el_rsc(:,i) - x(3)*v_el
            % or proj_el(:,j) = v1 + v13*x(1) + v12*x(2)

			v1_pel = v1-proj_el(:,j); no(1) = norm(v1_pel);
			v2_pel = v2-proj_el(:,j); no(2) = norm(v2_pel);
			v3_pel = v3-proj_el(:,j); no(3) = norm(v3_pel);
            if all(no>1e-7) % This is in case one electrode falls exactly on the vertex
                v1_pel = v1_pel/no(1);
                v2_pel = v2_pel/no(2);
                v3_pel = v3_pel/no(3);
				sum_angle = real(acos(v1_pel'*v2_pel) + ...
						 acos(v2_pel'*v3_pel) + ...
						 acos(v3_pel'*v1_pel)) ;
				if abs(sum_angle-2*pi)<1e-7
					is_in(j)=1 ;
				end
            else
                is_in(j)=1;
			end
        end
		which_tri = find(is_in==1) ;
		% if the projected electr. location is exactly on an edge, it is in 2 triangles
		%	=> keep the 1st one

		electrodes.vert(ii) = ind_v ;
		electrodes.tri(ii) = list_t(which_tri(1)) ;
        if flags.br_only
            % need to add the br_sc width here
            pos_el_mm(:,ii) = proj_el(:,which_tri(1))+v_el*flags.scbr_w;
        else
    		pos_el_mm(:,ii) = proj_el(:,which_tri(1)) ;
        end
	end

	electrodes.XYZmm = pos_el_mm;
    electrodes.el_loc = el_loc;
	electrodes.nr	 = Nel ;
    if flags.br_only
        electrodes.info	= ['Electr loc projected on brain surface + ',num2str(flags.scbr_w),' mm'];
    else
    	electrodes.info	= 'Electrode location projected on scalp surface' ;
    end
    electrodes.names = el_name';

	varargout{1} = electrodes;
    varargout{2} = flags;
%________________________________________________________________________
case 'genmesh'
%------------------------------------------------------------------------
% FORMAT head = spm_eeg_inv_model('GenMesh',Pvol,flags);
%
% Generate the tesselated surface of the 1 to 3 head volumes
% (brain-skull-scalp) from the binarized volumes.
% Input :
%   Pvol    : filenames of 3 (or 1-2) head volumes 'brain', 'skull' & 'scalp'.
%   flags    : various flags
%      * n       : (1x3) provides the number of vertices on each surface
%                   Npt = n(i)^2*5/4+2
%      * br_only : onlyt use the binarised brain volume (1) or scalp (0, default)
%      * q_elastm: Correct mesh using an elastic model (1, default), or not (0).
%      * q_meased: Measure edges of mesh (1), or not (0, default)
%      * q_4thpt : Determine 4th point central of each traingle (1), or not (0, default)
% Output :
%   head    : head structures with 1 to 3 tesselated surfaces
%------------------------------------------------------------------------
    fprintf('\n\n'), fprintf('%c','='*ones(1,80)), fprintf('\n')
    fprintf(['\tGenerate surface meshes from binarized volumes\n']);
    Pvol = varargin{1};
    
    def_flags = struct('Nvol',1,'n',[32 32 32],'br_only',0, ...
                        'q_elast_m',1,'q_meased',0,'q_4thpt',0);
    if nargin<3
        flags = def_flags;
    else
        flags = varargin{2};
        fnms  = fieldnames(def_flags);
		for i=1:length(fnms),
			if ~isfield(flags,fnms{i}),
                flags = setfield(flags,fnms{i},getfield(def_flags,fnms{i}));
            end
		end
    end

    [Centre_vx,Centre_mm]=spm_eeg_inv_model('CtrBin',Pvol(1,:));
    Nbin = size(Pvol,1);
    if length(flags.n)<Nbin
        n = flags.n(1)*ones(1,Nbin);
    else
        n = flags.n;
    end
    mesh_labels = strvcat('Tess. outer scalp','Tess. outer skull','Tess. inner skull');
    if flags.br_only
        head = spm_eeg_inv_model('TesBin',n(1),Centre_mm,Pvol,mesh_labels(3,:));
    elseif Nbin==2 % Brain & scalp only, no skull.
        head(1) = spm_eeg_inv_model('TesBin',n(1),Centre_mm,Pvol(1,:),mesh_labels(3,:));    
        head(2) = spm_eeg_inv_model('TesBin',n(2),Centre_mm,Pvol(2,:),mesh_labels(1,:));    
    else
        for ii=Nbin:-1:1
            head(ii) = spm_eeg_inv_model('TesBin',n(ii),Centre_mm,Pvol(ii,:),mesh_labels(ii,:));
        end
    end
    if flags.q_elast_m
    	fprintf(['\t\tCorrect position of vertices\n']);
        for ii=Nbin:-1:1
            head(ii) = spm_eeg_inv_model('ElastM',head(ii),Pvol(ii,:));
        end
    end
    if flags.q_meased
    	fprintf(['\t\tMeasure the edges length\n']);
        head(1).edges = struct('ed',[],'nr',[]);
        for ii=Nbin:-1:1
            head(ii) = spm_eeg_inv_model('MeasEd',head(ii));
        end
    end
    if flags.q_4thpt
		fprintf(['\t\tDetermine triangle''s "4th point"\n']);
        head(1).pt4 = [];
        for ii=Nbin:-1:1
            head(ii) = spm_eeg_inv_model('Tr4thPt',head(ii),Pvol(ii,:));
        end
    end
    varargout{1} = head;
    varargout{2} = Centre_mm;
    varargout{3} = flags;
    fprintf('%c','='*ones(1,80)), fprintf('\n')

%________________________________________________________________________
case 'tesbin'
%------------------------------------------------------------------------
% FORMAT ts = spm_eeg_inv_model('TesBin',n,ctr_vol,P,info);
%
% Generate a mesh covering a binarized volume
%   1. Generate a simple spherical mesh
%   2. The spherical mesh is projected radially on the bin volume
% Afterwards, "elastic" mesh correction can thus be useful to correct 
% some overlong edges.
%
% Input : 
%   n:      number of vertices on each surface Npt = n^2*5/4+2
%   ctr_vol:centre of bin volume for the radial projection (mm)
%   P:      filename of bin volume
%   info:   information string
%------------------------------------------------------------------------
	if nargin==1
        n = 32;
        P = spm_select(1,'outer*.img','Bin image to tessalate');
		[ctr_vol_vx,ctr_vol] = spm_eeg_inv_model('CtrBin',P) ;
		info = 'realistic model' ;
	elseif nargin==2
        n = varargin{1};
        P = spm_select(1,'outer*.img','Bin image to tessalate');
		[ctr_vol_vx,ctr_vol] = spm_eeg_inv_model('CtrBin',P) ;
		info = 'realistic model' ;
	elseif nargin==3
        n = varargin{1};
        ctr_vol = varargin{2};
        P = spm_select(1,'outer.img','Bin image to tessalate');
		info = 'realistic model' ;
	elseif nargin==4
        n = varargin{1};
        ctr_vol = varargin{2};
        P = varargin{3};
		info = 'realistic model' ;
	elseif nargin==5
        n = varargin{1};
        ctr_vol = varargin{2};
        P = varargin{3};
		info = varargin{4} ;
	elseif nargin>5
		error('Wrong input arguments for ''TesBin''.') ;
    end
    
	% Load volume information
	%------------------------
	Vv    = spm_vol(P);
    VOX   = sqrt(sum(Vv.mat(1:3,1:3).^2));
    
	% A few defintions
	%-----------------
	ho = 1 ; % Trilinear interpolation shoud be enough
	n_div = 1000; % # of sample in the radial direction
	trsh_vol = .9; % Thershold for surface detection
	rad = round(4/3*max(Vv.dim(1:3).*VOX/2)); % Radius (mm) of original sphere
	dr = rad/n_div; % size of sampling step (mm)
	d_li = 0:dr:rad; % Sampling location on radius (mm)
	nd_li = length(d_li); unit = ones(1,nd_li);
    
	% Create a tessalated sphere
	%---------------------------
	tsph = spm_eeg_inv_model('TesSph',n,rad);    % tesselated sphere of radius rad,
    vert = [tsph.vert ; ones(1,tsph.nr(1))] ;    % centered around [0 0 0]
	vert = spm_matrix([ctr_vol 0 pi/2 0])*vert ;        
	vert = vert(1:3,:) ;                                
    % Rotate the sphere by 90deg around y axis and centre sphere around the "centre" of the brain, 
    % All this in mm. Leave them as a 3xNvert matrix.
	
	srf_vert = zeros(3,tsph.nr(1)) ; % vertices at the surface of brain, in vx !
    spm_progress_bar('Init',tsph.nr(1),	['Generate ',info ],'Vertices projected');

	for i = 1:tsph.nr(1)
		or = ctr_vol'-vert(:,i) ; or = or/norm(or) ; % direction from the point toward the centre
		line = vert(:,i)*unit + or*d_li;
        line_vx = spm_eeg_inv_model('mm2vx',line,Vv.mat);
		val_line = spm_sample_vol(Vv,line_vx(1,:),line_vx(2,:),line_vx(3,:),ho) ;
		srf_vert(:,i) = line_vx(:,min(find(val_line>trsh_vol))) ;
                    % first point to intercept the surface
    	spm_progress_bar('Set',i);
    end
    spm_progress_bar('Clear');

	srf_vert_mm      = spm_eeg_inv_model('vx2mm',srf_vert,Vv.mat);
	surf.XYZmm       = srf_vert_mm(1:3,:);
	surf.tri         = tsph.tri;
	surf.nr          = tsph.nr;
	surf.M           = Vv.mat;
	surf.info.str    = info;
    surf.info.bin_fn = P;
    surf.Centre      = ctr_vol;
	varargout{1}     = surf;

%________________________________________________________________________
case 'elastm'
%------------------------------------------------------------------------
% FORMAT ts = spm_eeg_inv_model('ElastM',ts);
%
% Modify the the mesh in order to reduce overlong edges.
% The procedure uses an elastic model :
% At each vertex, the neighbouring triangles and vertices connected directly
% are considered.
% Each edge is considered elastic and can be lengthened or shortened,
% depending on their legnth.
% Refs: G.Taubin, A signal processing approach to fair surface design, 1995
% This is a non-shrinking smoothing algo.
%
% Input : 
%   ts:     tesselated surface
% Output :
%   ts:     tesselated surface with corrected mesh
%------------------------------------------------------------------------
    ts = varargin{1};

	M_con = sparse([ts.tri(1,:)';ts.tri(1,:)';ts.tri(2,:)';ts.tri(3,:)';ts.tri(2,:)';ts.tri(3,:)'], ...
                [ts.tri(2,:)';ts.tri(3,:)';ts.tri(1,:)';ts.tri(1,:)';ts.tri(3,:)';ts.tri(2,:)'], ...
                ones(ts.nr(2)*6,1),ts.nr(1),ts.nr(1));
            % Connection vertex-to-vertex
                      
	kpb = .1; % Cutt-off frequency
	lam = .5; mu = lam/(lam*kpb-1); % Parameters for elasticity.
	N = 25; % Number of smoothing steps, the larger, the smoother
	
	XYZmm = ts.XYZmm;
	for j=1:N
        XYZmm_o = zeros(3,ts.nr(1)) ;
		XYZmm_o2 = zeros(3,ts.nr(1)) ;
		for i=1:ts.nr(1)
		    ln = find(M_con(:,i));
            d_i = sqrt(sum((XYZmm(:,ln)-XYZmm(:,i)*ones(1,length(ln))).^2));
            w_i = d_i/sum(d_i);
            XYZmm_o(:,i) = XYZmm(:,i) + ...
                lam * sum((XYZmm(:,ln)-XYZmm(:,i)*ones(1,length(ln))).*(ones(3,1)*w_i),2);
		end
		for i=1:ts.nr(1)
		    ln = find(M_con(:,i));
            d_i = sqrt(sum((XYZmm(:,ln)-XYZmm(:,i)*ones(1,length(ln))).^2));
            w_i = d_i/sum(d_i);
            XYZmm_o2(:,i) = XYZmm_o(:,i) + ...
                mu * sum((XYZmm_o(:,ln)-XYZmm_o(:,i)*ones(1,length(ln))).*(ones(3,1)*w_i),2);
		end
        XYZmm = XYZmm_o2;
	end
	varargout{1} = ts;
	varargout{1}.XYZmm = XYZmm;
%________________________________________________________________________
case 'meased'
%------------------------------------------------------------------------
% FORMAT ts = spm_eeg_inv_model('MeasEd',ts);
%
% This measures the edges of the mesh.
%   Input : ts
%   Output : ts with measured edges
%
% Format of the '.edges.' sub-structure :
% 	ts.edges.ed :
% 		ed(:,1) = length,
%		ed(:,[2 3]) = indexes of 2 points,
%		ed(:,[4 5]) = indexes of 2 triangles.
%	ts.edges.nr :
%		nr(1) : number of edges,
%		nr([2 3]) : mean and std of edges length.
%------------------------------------------------------------------------

ts = varargin{1}; XYZmm = ts.XYZmm;
Ned = 3*ts.nr(2)/2; edges = zeros(Ned,5); i_ed = 1 ;
for i=1:ts.nr(2) % Go through all the triangles
	e1 = ts.tri(1:2,i); % 3 edges of the ith triangle (counter-clock wise)
	e2 = ts.tri(2:3,i);
	e3 = ts.tri([3 1],i);
	p1 = find((edges(:,2)==e1(2)) && (edges(:,3)==e1(1))); % check if edge was already seen.
	p2 = find((edges(:,2)==e2(2)) && (edges(:,3)==e2(1)));
	p3 = find((edges(:,2)==e3(2)) && (edges(:,3)==e3(1)));
	if isempty(p1) % new edge
		edges(i_ed,:) = [norm(XYZmm(:,e1(1))-XYZmm(:,e1(2))) e1' i 0] ;
		i_ed = i_ed+1;
	else % 2nd triangle for this edge
		edges(p1,5) = i;
	end
	if isempty(p2)
		edges(i_ed,:) = [norm(XYZmm(:,e2(1))-XYZmm(:,e2(2))) e2' i 0] ;
		i_ed = i_ed+1 ;
	else
		edges(p2,5) = i;
	end
	if isempty(p3)
		edges(i_ed,:) = [norm(XYZmm(:,e3(1))-XYZmm(:,e3(2))) e3' i 0] ;
		i_ed = i_ed+1 ;
	else
		edges(p3,5) = i;
	end
end
ts.edges.ed = edges; ts.edges.nr = [Ned mean(edges(:,1)) std(edges(:,1))];
varargout{1} = ts;

%________________________________________________________________________
case 'tr4thpt'
%------------------------------------------------------------------------
% FORMAT ts = spm_eeg_inv_model('Tr4thPt',ts);
%
% Finds the 4th' point of each triangle:
% a flat triangle is only an approximation of the real underlying volume surface, 
% a 4th point taken at the centre of the triangle is found by approximating
% the neighbouring triangles in a sphere, to provide an idea of the local curvature.
% This is used in the BEM solution of the forward problem, to sort out 
% the problem of the auto-solid angle.
%
% Input : 
%   ts:     tesselated surface
% Output :
%   ts:     tesselated surface with 4th point at each triangle
%------------------------------------------------------------------------
	if nargin==1
        error('Wrong input. You should pass a tess. surface')
	end
    ts = varargin{1};
	T_con = sparse([ts.edges.ed(:,4),ts.edges.ed(:,5)],[ts.edges.ed(:,5),ts.edges.ed(:,4)], ...
                ones(ts.edges.nr(1)*2,1),ts.nr(2),ts.nr(2),ts.edges.nr(1)*2);
            % Connection Triangle-to-triangle
	
	XYZmm = ts.XYZmm;
	
	ts.pt4.XYZmm = zeros(3,ts.nr(2)) ;
	for i=1:ts.nr(2)
        l_tr = find(T_con(:,i));
        l_vt = ts.tri(:,l_tr); l_vt = sort(l_vt(:));
        l_vt = [l_vt(find(diff(l_vt)));l_vt(end)]; % Remov vertices listed twice.
        vi = XYZmm(:,l_vt); A = [2*vi' -ones(6,1)]; b = sum(vi.^2)';
        x = A\b;
        vc = x(1:3);
        R = sqrt(norm(vc)^2-x(4));
        vct = sum(XYZmm(:,ts.tri(:,i)),2)/3; or = vct-vc; or = or/norm(or);
        ts.pt4.XYZmm(:,i) = vc+R*or;
    end
    varargout{1} = ts;
    
%________________________________________________________________________
case 'tessph'    
%------------------------------------------------------------------------
% FORMAT tsph = spm_eeg_inv_model('TesSph',n,r);
%
% Generate a structure 'tsph' containing a tesselated sphere.
% Input : 
%   n : number of 'latitude' divisions on the sphere. It MUST be even!
%   r : radius of the sphere
% Output : 
%   tsph .vert : vertices coordinates (3 x Nvert)
%		 .tri  : triangle patches (3 x Ntri)
%		 .info : info string
% 		8  -> 82   points -> 160 triangles
% 		10 -> 127  points -> 250 triangles
% 		12 -> 182  points -> 360 triangles
% 		14 -> 247  points -> 490 triangles
% 		16 -> 322  points -> 640 triangles
% 		18 -> 407  points -> 810 triangles
% 		20 -> 502  points -> 1000 triangles
% 		22 -> 607  points -> 1210 triangles
% 		24 -> 722  points -> 1440 triangles
% 		26 -> 847  points -> 1690 triangles
% 		28 -> 982  points -> 1960 triangles
% 		30 -> 1127 points -> 2250 triangles
% 		32 -> 1282 points -> 2560 triangles
% 		34 -> 1447 points -> 2890 triangles
% 		36 -> 1622 points -> 3240 triangles
% 		40 -> 2002 points -> 4000 triangles
% 		42 -> 2207 points -> 4410 triangles
% 		44 -> 2422 points -> 4840 triangles
% 		46 -> 2647 points -> 5290 triangles
% 		48 -> 2882 points -> 5760 triangles
% 		54 -> 3647 points -> 7290 triangles
% 		56 -> 3922 points -> 7840 triangles
% 		58 -> 4207 points -> 8410 triangles
% 		60 -> 4502 points -> 9000 triangles
% 		62 -> 4807 points -> 9610 triangles
% 		64 -> 5122 points -> 10240 triangles.
%------------------------------------------------------------------------
    if nargin==1 
		n = 20 ;
		r = 1 ;
		disp(['default values : n=',num2str(n),' , r=',num2str(r)]) ;
	elseif nargin==2
            n = varargin{1}
            r = 1 ;
		    disp(['default value :  r=',num2str(r)]) ;
	elseif nargin==3
            n = varargin{1};
            r = varargin{2};
	else
        error('Wrong input format 2')
	end
	
	[X_sp,Y_sp,Z_sp,npt_tr,npt_qr,npt_sph] = point(n,r) ;
	[ind_tr] = tri(npt_tr,npt_qr,npt_sph) ;
	
	tsph.vert = [X_sp ; Y_sp ; Z_sp] ;
	tsph.tri = ind_tr ;
	tsph.nr = [length(X_sp) size(ind_tr,2)] ;
	tsph.info = ['Tes sph : n=',num2str(n),', r=',num2str(r)] ;
	varargout{1} = tsph;
%__________________________________________________________________________
case 'erodegrow'    
%------------------------------------------------------------------------
%
% FORMAT [Pout/val] = spm_eeg_inv_model('ErodeGrow',P/val,ne,ng,thr_im)
% 
% It erodes then grows an image after thresholding it.
% Inputs : P/val,ne,ng,thr_im
%     P : file name of image to erode/grow
%     val : full volume of image (as loaded by spm_read_vols)
%     ne : nr of steps for erosion, default 3
%     ng : nr of steps for growing, default 6
%     thr_im : threshold value to apply, default .8
%                (compared to the maximal value of image)
% Output : Pout/val
%     Pout : name of the file generated
%     val : full value of eroded-grown image
%     
% Notes:
% - if 2 file names are specified in P[input], the 2nd one is used as the name
%   of the generated file otherwise the new filename is creted from P as
%           [P,'_e',num2str(ne),'g',num2str(ng),'.img']
% - If a file name is passed the output is a filename.
%   If a matrix of values is passed, the output is a matrix of values.
%------------------------------------------------------------------------
fl_rvol = 0; % Need to load (1) or not (0) the volume from a file
if nargin<2
    P = spm_select(1,'image','Image to erode-grow');
    nPin = size(P,1);
    Vp = spm_vol(P);
    fl_rvol = 1;
    dim = V.dim;
    cr_file = 1; % Create (1) or not (0) a file at the end of the process.
end

if length(varargin)>=1
    if ischar(varargin{1})
        P = varargin{1};
        nPin = size(P,1);
        Vp = spm_vol(P(1,:));
        fl_rvol = 1;
        dim = Vp.dim;
        cr_file = 1;
    elseif isa(varargin{1},'uint8')
        val = varargin{1};
        dim = size(val);
        cr_file = 0;
    elseif isnumeric(varargin{1})
        val = uint8(varargin{1});
        dim = size(val);
        cr_file = 0;
    else
        error('Wrong data format');
    end
end

if fl_rvol
    val = spm_loaduint8(Vp);
end

ne = 1; ng = 1; thr_im = .1 ;
if length(varargin)==4
    ne =  varargin{2}; 
    ng =  varargin{3}; 
    thr_im =  varargin{4};
else
    disp('Using default erode/grow parameters: ne=1, ng=1, thr=.1 !')
end
p1 = uint8(zeros(1,dim(2),dim(3)));
p2 = uint8(zeros(dim(1),1,dim(3)));
p3 = uint8(zeros(dim(1),dim(2),1));

p1_3 = uint8(zeros(1,dim(2),dim(3)-1));
p2_1 = uint8(zeros(dim(1)-1,1,dim(3)));
p3_2 = uint8(zeros(dim(1),dim(2)-1,1));

% Erosion !
% val = val>thr_im*max(max(max(abs(val)))) ;
val = uint8(val>thr_im*double(max(val(:)))) ; 
    % added uint8 otherwise is considered as a logical variable 
    %   -> no math operation (addition) possible
Nbin = length(find(val(:)));
fprintf('\tOriginal number of voxels %d .\n',Nbin)
for i=1:ne
    temp = val ...
        + cat(1,p1,val(1:end-1,:,:)) + cat(1,val(2:end,:,:),p1) ...
        + cat(2,p2,val(:,1:end-1,:)) + cat(2,val(:,2:end,:),p2) ...
        + cat(3,p3,val(:,:,1:end-1)) + cat(3,val(:,:,2:end),p3) ...
        + cat(1,p1,cat(2,p2_1,val(1:end-1,1:end-1,:))) ...
        + cat(1,p1,cat(2,val(1:end-1,2:end,:),p2_1)) ...
        + cat(1,cat(2,p2_1,val(2:end,1:end-1,:)),p1) ...
        + cat(1,cat(2,val(2:end,2:end,:),p2_1),p1) ...        
        + cat(2,p2,cat(3,p3_2,val(:,1:end-1,1:end-1))) ...
        + cat(2,p2,cat(3,val(:,1:end-1,2:end),p3_2)) ...
        + cat(2,cat(3,p3_2,val(:,1:end-1,1:end-1)),p2) ...
        + cat(2,cat(3,val(:,1:end-1,2:end),p3_2),p2) ...        
        + cat(3,p3,cat(1,p1_3,val(1:end-1,:,1:end-1))) ...
        + cat(3,p3,cat(1,val(2:end,:,1:end-1),p1_3)) ...
        + cat(3,cat(1,p1_3,val(1:end-1,:,1:end-1)),p3) ...
        + cat(3,cat(1,val(2:end,:,1:end-1),p1_3),p3) ;
    val = uint8(temp>17) ;
 %   val = temp>6 ;
    Nbin = length(find(val(:)));
    fprintf('\tErosion step %d , %d voxels left.\n',i,Nbin)
end

%Ve = Vp;
%Ve.fname = [spm_str_manip(Vp.fname,'r'),'_e',num2str(ne),'.img']
%Vout_e = spm_write_vol(Ve,val)

% Growing bit !
for i=1:ng
    temp = val ...
        + cat(1,p1,val(1:end-1,:,:)) + cat(1,val(2:end,:,:),p1) ...
        + cat(2,p2,val(:,1:end-1,:)) + cat(2,val(:,2:end,:),p2) ...
        + cat(3,p3,val(:,:,1:end-1)) + cat(3,val(:,:,2:end),p3) ...
        + cat(1,p1,cat(2,p2_1,val(1:end-1,1:end-1,:))) ...
        + cat(1,p1,cat(2,val(1:end-1,2:end,:),p2_1)) ...
        + cat(1,cat(2,p2_1,val(2:end,1:end-1,:)),p1) ...
        + cat(1,cat(2,val(2:end,2:end,:),p2_1),p1) ...        
        + cat(2,p2,cat(3,p3_2,val(:,1:end-1,1:end-1))) ...
        + cat(2,p2,cat(3,val(:,1:end-1,2:end),p3_2)) ...
        + cat(2,cat(3,p3_2,val(:,1:end-1,1:end-1)),p2) ...
        + cat(2,cat(3,val(:,1:end-1,2:end),p3_2),p2) ...        
        + cat(3,p3,cat(1,p1_3,val(1:end-1,:,1:end-1))) ...
        + cat(3,p3,cat(1,val(2:end,:,1:end-1),p1_3)) ...
        + cat(3,cat(1,p1_3,val(1:end-1,:,1:end-1)),p3) ...
        + cat(3,cat(1,val(2:end,:,1:end-1),p1_3),p3) ;
    val = uint8(temp>1) ;
    Nbin = length(find(val(:)));
    fprintf('\tGrowing step %d , %d voxels left.\n',i,Nbin)
end

if cr_file
    Vout = Vp;
%     val = val.*uint8(255);
    Vout.pinfo = [1/255 0 0]';
    if nPin==1
        Vout.fname = [spm_str_manip(Vp.fname,'r'),'_e',num2str(ne),'g',num2str(ng),'.img'];
    else
        Vout.fname = P(2,:);
    end
    spm_write_vol(Vout,val);
    varargout{1} = Vout.fname ;
else
    varargout{1} = val ;
end



%________________________________________________________________________
otherwise,
	warning('Unknown action string')
end;

return;

%________________________________________________________________________
%________________________________________________________________________
%________________________________________________________________________
%________________________________________________________________________
%
% SUBFUNCTIONS
%________________________________________________________________________
%________________________________________________________________________
%________________________________________________________________________
% FORMAT [X_sp,Y_sp,Z_sp,npt_tr,npt_qr,npt_sph] = point(n,r)
% Generate a set of points on a sphere.
% n = the number of latitudes on the sphere, it MUST be even!
% r = radius of the sphere, 1 by default.
%------------------------------------------------------------------------
function [X_sp,Y_sp,Z_sp,npt_tr,npt_qr,npt_sph] = point(n,r)
d_r = pi/180; dth = 180/n*d_r ;
X_qr=[] ; Y_qr=[] ; Z_qr=[] ; X_sp=[] ; Y_sp=[] ; Z_sp=[] ;

% Coord of pts on a 5th of sphere WITHOUT the top and bottom points
for i=1:n/2         % Upper part
	the = i*dth ;
	cth = cos(the) ; sth = sin(the) ;
	dphi = 72/i*d_r ;
	for j=1:i
		phi = (j-1)*dphi ;
		cph = cos(phi) ; sph = sin(phi) ;
		X_qr = [X_qr sth*cph] ;
		Y_qr = [Y_qr sth*sph] ;
		Z_qr = [Z_qr cth] ;
	end
end
for i=(n/2+1):(n-1) % Lower part
	the = i*dth ;
	cth = cos(the) ; sth = sin(the) ;
	ii = n-i ;
	dphi = 72/ii*d_r ;
	for j=1:ii
		phi = (j-1)*dphi ;
		cph = cos(phi) ; sph = sin(phi) ;
		X_qr = [X_qr sth*cph] ;
		Y_qr = [Y_qr sth*sph] ;
		Z_qr = [Z_qr cth] ;
	end
end
% Each part is copied with a 72deg shift
dphi=-72*d_r ;
for i=0:4
	phi=i*dphi ;
	cph=cos(phi) ;
	sph=sin(phi) ;
	X_sp=[X_sp cph*X_qr+sph*Y_qr] ;
	Y_sp=[Y_sp -sph*X_qr+cph*Y_qr] ;
	Z_sp=[Z_sp Z_qr] ;
end
% add top and bottom points
X_sp=[0 X_sp 0]*r ;
Y_sp=[0 Y_sp 0]*r ;
Z_sp=[1 Z_sp -1]*r ;

npt_tr = [[0:n/2] [n/2-1:-1:0]]; % Nbr of points on each slice for a 5th of sph.
npt_qr = sum(npt_tr) ; % total nbr of pts per 5th of sphere
npt_sph=5/4*n^2+2 ; % or 5*npt_qr+2 ; % total nbr of pts on sph
return

%________________________________________________________________________
% FORMAT [ind_tr] = tri(npt_tr,npt_qr,npt_sph)
%------------------------------------------------------------------------
% Order the pts created by the function 'point' into triangles
% Each triangle is represented by 3 indices correponding to the XYZ_sp
%
% This works only because I know in wich order the vertices are generated
% in the 'point' function.
%------------------------------------------------------------------------
function [ind_tr] = tri(npt_tr,npt_qr,npt_sph)
n = sqrt(4/5*(npt_sph-2)) ;

% First 5th of sphere only
for i=0:(n/2-1)         % upper half
	if i==0		    	% 1st triangle at the top
		ind_qr=[1 2 2+npt_qr]' ;
	else		    	% other triangles
		for j=0:npt_tr(i+1)
			if j==npt_tr(i+1)-1			% last but 1 pt on slice
				x1=sum(npt_tr(1:i))+j+2 ;
				x12=x1+npt_tr(i+1) ;
				x13=x12+1 ;
				x22=x13 ;
				x23=sum(npt_tr(1:i))+2+npt_qr ;
				ind_qr=[ind_qr [x1 x12 x13]' [x1 x22 x23]'] ;
			elseif j==npt_tr(i+1)		% last pt on slice
				x1=sum(npt_tr(1:i))+2+npt_qr ;
				x2=sum(npt_tr(1:(i+1)))+j+2 ;
				x3=x1+npt_tr(i+1) ;
				ind_qr=[ind_qr [x1 x2 x3]'] ;
			else        			% other pts on slice
				x1=sum(npt_tr(1:i))+j+2 ;
				x12=x1+npt_tr(i+1) ;
				x13=x12+1 ;
				x22=x13 ;
				x23=x1+1 ;
				ind_qr=[ind_qr [x1 x12 x13]' [x1 x22 x23]'] ;
			end
		end
	end
end
for i=(n/2+1):n         % lower half
	if i==n	        	% last triangle at the bottom
		ind_qr=[ind_qr [5*npt_qr+2 2*npt_qr+1 npt_qr+1]' ] ;
	else            	% other triangles
		for j=0:npt_tr(i+1)
			if j==npt_tr(i+1)-1			% last but 1 pt on slice
				x1=sum(npt_tr(1:i))+j+2 ;
				x12=x1-npt_tr(i) ;
				x13=x12+1 ;
				x22=x13 ;
				x23=sum(npt_tr(1:i))+2+npt_qr ;
				ind_qr=[ind_qr [x1 x13 x12]' [x1 x23 x22]'] ;
			elseif j==npt_tr(i+1)		% last pt on slice
				x1=sum(npt_tr(1:i))+2+npt_qr ;
				x2=sum(npt_tr(1:(i-1)))+j+2 ;
				x3=x1-npt_tr(i) ;
				ind_qr=[ind_qr [x1 x3 x2]'] ;
			else            			% other pts on slice
				x1=sum(npt_tr(1:i))+j+2 ;
				x12=x1-npt_tr(i) ;
				x13=x12+1 ;
				x22=x13 ;
				x23=x1+1 ;
				ind_qr=[ind_qr [x1 x13 x12]' [x1 x23 x22]'] ;
			end
		end
	end
end


ntr_qr=size(ind_qr,2) ; % nbr of triangle per 5th of sphere
[S_i,S_j]=find(ind_qr==1) ;
[B_i,B_j]=find(ind_qr==npt_sph) ;
[qs_i,qs_j]=find((ind_qr>(npt_qr+1))&(ind_qr<(3*npt_qr))) ;
ind_tr=[] ;

% shift all indices to cover the sphere
for i=0:4
	ind_tr = [ind_tr (ind_qr+i*npt_qr)] ;
	ind_tr(S_i,S_j+i*ntr_qr) = 1 ;
	ind_tr(B_i,B_j+i*ntr_qr) = npt_sph ;
end
for i=1:size(qs_i,1)
	ind_tr(qs_i(i),qs_j(i)+4*ntr_qr) = ...
		ind_tr(qs_i(i),(qs_j(i)+4*ntr_qr))-5*npt_qr ;
end
ntr_sph = size(ind_tr,2) ;% nbr of triangles on the sphere
	% or = 5/2*n^2
return

% %________________________________________________________________________
% % FORMAT [udat] = loaduint8(V)
% %------------------------------------------------------------------------
% % Load data from file indicated by V into an array of unsigned bytes.
% %------------------------------------------------------------------------
% 
% function udat = loaduint8(V)
% % Load data from file indicated by V into an array of unsigned bytes.
% 
% if size(V.pinfo,2)==1 && V.pinfo(1) == 2,
% 	mx = 255*V.pinfo(1) + V.pinfo(2);
% 	mn = V.pinfo(2);
% else,
% 	spm_progress_bar('Init',V.dim(3),...
% 		['Computing max/min of ' spm_str_manip(V.fname,'t')],...
% 		'Planes complete');
% 	mx = -Inf; mn =  Inf;
% 	for p=1:V.dim(3),
% 		img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
% 		mx  = max([max(img(:))+paccuracy(V,p) mx]);
% 		mn  = min([min(img(:)) mn]);
% 		spm_progress_bar('Set',p);
% 	end;
% end;
% spm_progress_bar('Init',V.dim(3),...
% 	['Loading ' spm_str_manip(V.fname,'t')],...
% 	'Planes loaded');
% 
% udat = uint8(0);
% udat(V.dim(1),V.dim(2),V.dim(3))=0;
% rand('state',100);
% for p=1:V.dim(3),
% 	img = spm_slice_vol(V,spm_matrix([0 0 p]),V.dim(1:2),1);
% 	acc = paccuracy(V,p);
% 	if acc==0,
% 		udat(:,:,p) = uint8(round((img-mn)*(255/(mx-mn))));
% 	else,
% 		% Add random numbers before rounding to reduce aliasing artifact
% 		r = rand(size(img))*acc;
% 		udat(:,:,p) = uint8(round((img+r-mn)*(255/(mx-mn))));
% 	end;
% 	spm_progress_bar('Set',p);
% end;
% spm_progress_bar('Clear');
% return;
% 
% function acc = paccuracy(V,p)
% % if ~spm_type(V.dim(4),'intt'),
% if ~spm_type(V.dt(1),'intt'),
% 	acc = 0;
% else,
% 	if size(V.pinfo,2)==1,
% 		acc = abs(V.pinfo(1,1));
% 	else,
% 		acc = abs(V.pinfo(1,p));
% 	end;
% end;
%________________________________________________________________________
% FORMAT savefields(fnam,p)
%------------------------------------------------------------------------
% save the fields of a structure 'p' in a file 'fnam'
%------------------------------------------------------------------------
function savefields(fnam,p)
if length(p)>1, error('Can''t save fields.'); end;
fn = fieldnames(p);
for i=1:length(fn),
    eval([fn{i} '= p.' fn{i} ';']);
end;
if spm_matlab_version_chk('7') >= 0,
    save(fnam,'-V6',fn{:});
else
    save(fnam,fn{:});
end;
return;
