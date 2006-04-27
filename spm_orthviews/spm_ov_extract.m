function ret = spm_ov_extract(varargin)
% Wrapper for tbxvol_extract
% This function provides an interface to tbxvol_extract from
% an spm_orthviews display. The main features include:
% ROI selection alternatives
% - cross hair selection for start and end point of a profile line
% (after selection, the line will be added as a yellow blob to the
% displayed image)
% - ROI specification based on displayed blobs (1 ROI per blob set)
% - ROI specification from active ROI tool
% - ROI specification from saved files (1 ROI per file)
% - Current crosshair position
% - Sphere at current crosshair position
% Extraction alternatives
% - Raw and adjusted data from SPM analysis (SPM.mat)
% - Raw data from image series
% - Raw data from current image only
% Plotting alternatives
% - no plot (just save extracted data)
% - plot over voxels (1 line per image)
% - plot over images (1 line per voxel)
% - average over voxels/images (average line and scatter plot of
%   individual values)
% - if DTI toolbox installed: average over voxels and diffusion
%   weighting
% If data is extracted from the current image, then only plotting
% over voxels is available.
%
% This routine is a plugin to spm_orthviews for SPM5. For general help about
% spm_orthviews and plugins type
%             help spm_orthviews
% at the matlab prompt.
%_______________________________________________________________________
%
% @(#) $Id: spm_ov_extract.m,v 1.14 2006/02/07 09:00:10 glauche Exp $

rev = '$Revision: 1.14 $';

global st;
global defaults;
if isempty(st)
  error(['%s: This routine can only be called as a plugin for' ...
	 ' spm_orthviews!'], mfilename);
end;

if nargin < 2
  error(['%s: Wrong number of arguments. Usage:' ...
	 ' spm_orthviews(''extract'', cmd, volhandle, varargin)'], ...
	mfilename);
end;

cmd = lower(varargin{1});
volhandle = varargin{2};
files = [];

switch cmd
  
        %-------------------------------------------------------------------------
        % Context menu and callbacks
case 'context_menu'
        addpath(fullfile(spm('dir'),'toolbox','Volumes', ...
                         'Single_Volumes'));
        if ~any(exist('tbxvol_extract')==[2:6])
                warning([mfilename ':init'], 'function tbxvol_extract not found!');
                return;
        end;
        item0 = uimenu(varargin{3}, 'Label', 'Extract+Plot Data', ...
                       'Tag', sprintf('%s_0_%d', upper(mfilename),...
                                      num2str(volhandle)));
        item1 = uimenu(item0, 'Label','Time Series (SPM): - Raw and Fitted Data',... 
                       'Callback', sprintf(...
                           ['feval(''%s'',''context_init'', %d,'...
                            ' ''SPMboth'', ''ser'');'], ...
                           mfilename, volhandle),...
                       'Tag', sprintf('%s_0_%d', upper(mfilename),...
                                      num2str(volhandle)));
        item1 = uimenu(item0, 'Label','Image Series: - Raw Data',... 
                       'Callback', sprintf(...
                           ['feval(''%s'',''context_init'', %d,'...
                            ' ''ROIfiles'',''ser'');'], ...
                           mfilename, volhandle),...
                       'Tag', sprintf('%s_0_%d', upper(mfilename),...
                                      num2str(volhandle)));
        item1 = uimenu(item0, 'Label','Current Image: - Raw Data',... 
                       'Callback', sprintf(...
                           ['feval(''%s'',''context_init'', %d,'...
                            ' ''ROIfiles'',''cur'');'], ...
                           mfilename, volhandle),...
                       'Tag', sprintf('%s_0_%d', upper(mfilename),...
                                      num2str(volhandle)));
        item10 = uimenu(item0, 'Label', 'Help', 'Callback', ...
                        ['feval(''spm_help'',''' mfilename ''');']);
case 'context_init'
        mode = lower(varargin{3});
        possel = {'Sphere at Current Point','Current Point','ROI file(s)'};
        posselsrc = {'sphere','point','froi'};
        if isfield(st.vols{volhandle},'roi')
                possel{end+1} = 'Current ROI';
                posselsrc{end+1} = 'croi';
        elseif isfield(st.vols{volhandle},'blobs')
                possel{end+1} = 'Current blobs';
                posselsrc{end+1} = 'blobs';
        end;
        if strcmp(mode,'roifiles')
                possel{end+1} = 'Trajectory through image';
                posselsrc{end+1} = 'line';
        end;
        if strcmp(varargin{4},'cur')
                bch.src.srcimgs = {st.vols{volhandle}.fname};
                plotsel = {'Don''t plot', 'Plot over voxels - no averaging'};
                plotselsrc = [0 1];
        else
                if strcmp(mode,'roifiles')
                        bch.src.srcimgs = cellstr(spm_select([1 inf],'image', ...
                                                             'Images to Extract Data from'));
                else
                        bch.src.srcspm = {spm_select([1 1],'.*SPM.*\.mat', ...
                                                     'SPM Analysis to Extract Data from')};
                end;
                plotsel = {'Don''t plot', 'Plot over voxels - no averaging',...
                           'Plot over images - no averaging',...
                           'Plot over voxels - Average over images'...
                           'Plot over images - Average over voxels'};
                plotselsrc = 0:3;
                if (exist(fullfile(spm('dir'),'toolbox',...
                                   'Diffusion','Helpers','dti_extract_dirs.m'))|...
                    exist('dti_extract_dirs'))
                        addpath(fullfile(spm('dir'),'toolbox','Diffusion','Helpers'))
                        plotsel{end+1} = 'DTI: Average over voxels and diffusion weighting';
                        plotselsrc(end+1)=4;
                end;
        end;
        possel = char(spm_input('Voxels to extract from', '+1', 'm', possel, ...
                                posselsrc));
        plotsel = spm_input('Plot data','!+1','m',...
                            plotsel, 0:(numel(plotsel)-1), 0);
        bch.interp = spm_input('Interpolation hold','+1','m',...
                               ['Nearest Neighbour|Trilinear|2nd Degree B-Spline|'...
                            '3rd Degree B-Spline|4th Degree B-Spline|5th Degree B-Spline|'...
                            '6th Degree B-Spline|7th Degree B-Spline'],...
                               [0 1 2 3 4 5 6 7],0);
        switch possel
        case 'croi'
                bch.roispec{1}.roilist = inv(st.vols{volhandle}.roi.Vroi.mat)*...
                    [st.vols{volhandle}.roi.xyz; ...
                     ones(1,size(st.vols{volhandle}.roi.xyz,2))];
                bch.roispec{1}.roilist = bch.roispec{1}.roilist(1:3,:);
        case 'froi'
                [froiname sts] = spm_select([1 Inf], 'image', 'Select ROI file(s)');
                if ~sts
                        warning('%s: No ROI file selected - exiting\n', mfilename);
                        return;
                end;
                for k = 1:size(froiname,1)
                        bch.roispec{k}.srcimg = {froiname(k,:)};
                end;
        case 'blobs'
                nroi = numel(st.vols{volhandle}.blobs);
                for k = 1:nroi
                        if isstruct(st.vols{volhandle}.blobs{k}.vol)
                                dat=spm_read_vols(st.vols{volhandle}.blobs{k}.vol);
                        else
                                dat=st.vols{volhandle}.blobs{k}.vol;
                        end;
                        dat(isnan(dat))=0;
                        [x y z] = ind2sub(size(dat), find(dat));
                        bch.roispec{k}.roilist = st.vols{volhandle}.blobs{k}.mat*...
                            [x(:)'; y(:)'; z(:)'; ones(size(x(:)'))];
                        bch.roispec{k}.roilist = bch.roispec{k}.roilist(1:3,:);
                end;
        case 'point',
                bch.roispec{1}.roilist = spm_orthviews('pos');
        case 'sphere',
                bch.roispec{1}.roisphere.roicent = spm_orthviews('pos');
                bch.roispec{1}.roisphere.roirad  = spm_input(...
                    'Radius of sphere (mm)', '!+1', 'e', '10', [1 1]);
        case 'line'
                msgbox(['Place cross hair at start position and press middle mouse' ...
                        ' button'],'Extract data');
                st.vols{volhandle}.extract=struct('cb',[], 'plotsel',plotsel,...
                                                  'bch',bch, 'mode',mode);
                clear bch;
                for k = 1:3
                        st.vols{volhandle}.extract.cb{k} = ...
                            get(st.vols{volhandle}.ax{k}.ax, 'ButtonDownFcn');
                        set(st.vols{volhandle}.ax{k}.ax,...
                            'ButtonDownFcn',...
                            ['switch get(gcf,''SelectionType'')',...
                             'case ''normal'', spm_orthviews(''Reposition'');',...
                             'case ''extend'', spm_orthviews(''extract'',''getline'',', ...
                             num2str(volhandle), ');',...
                             'case ''alt'', spm_orthviews(''context_menu'',''ts'',1);',...
                             'end;']);
                end;
        end;
case 'getline'
        if ~isfield(st.vols{volhandle}.extract,'start')
                st.vols{volhandle}.extract.start = spm_orthviews('pos');
                msgbox(['Place cross hair at end position and press middle mouse' ...
                        ' button'],'Extract data');
        else
                bch = st.vols{volhandle}.extract.bch;
                bch.roispec{1}.roiline.roistart = st.vols{volhandle}.extract.start;
                bch.roispec{1}.roiline.roiend = spm_orthviews('pos');
                bch.roispec{1}.roiline.roistep = spm_input('Sampling Distance (mm)', ...
                                                           '!+1', 'e', '1', [1 1]);
                for k = 1:3
                        set(st.vols{volhandle}.ax{k}.ax, 'ButtonDownFcn',...
                                          st.vols{volhandle}.extract.cb{k});
                end;
                mode = st.vols{volhandle}.extract.mode;
                plotsel = st.vols{volhandle}.extract.plotsel;
                st.vols{volhandle} = rmfield(st.vols{volhandle},'extract');
        end;
case 'redraw'
        % do nothing
otherwise    
        fprintf('spm_orthviews(''extract'', ...): Unknown action %s', cmd);
end;
if exist('bch','var')
        ext=tbxvol_extract(bch);
        assignin('base','ext',ext);
        disp('extracted data saved to workspace variable ''ext''');
        disp(ext);
        if strcmp(cmd,'getline')
                spm_orthviews('addcolouredblobs', volhandle, ...
                              inv(st.vols{volhandle}.mat)*...
                              [ext.posmm;ones(1,  size(ext.posmm,2))],...
                              ones(1,size(ext.posmm,2)), st.vols{volhandle}.mat,...
                              [1 1 0]);
                spm_orthviews('redraw');
        end;
        for sub = 1:size(ext,2)
                switch plotsel
                case 1,
                        dirs.dirsel='voxno';
                case 2,
                        dirs.dirsel='imgno';
                case 3,
                        dirs.dirsel='voxavg';
                case 4,
                        dirs.dirsel='imgavg';
                case 5,
                        dirsbch.saveinf = 0;
                        dirsbch.ref.refscanner = 1;
                        dirsbch.sep = 1;
                        dirsbch.dtol = ...
                            defaults.tools.vgtbx_Diffusion.extract_dirs.dtol;
                        dirsbch.ltol = ...
                            defaults.tools.vgtbx_Diffusion.extract_dirs.ltol;
                        try
                                load(bch.src.srcspm{1});
                                dirsbch.srcimgs=SPM.xY.P;
                                dirs = dti_extract_dirs(dirsbch);
                        catch
                                try 
                                        dirsbch.srcimgs = bch.src.srcimgs;
                                        dirs = dti_extract_dirs(dirsbch);
                                catch
                                        disp('Cannot extract DTI information for plotting.');
                                end;
                        end;
                end;
                if plotsel
                        if ~isempty(ext(1,sub).raw)
                                try 
                                        data = cat(1,ext(:,sub).raw);
                                catch
                                        fprintf(['Can not plot data from ' ...
                                                'images with different ' ...
                                                'voxel sizes and orientations.\n']);
                                        error(lasterr);
                                end;
                                extract_fig(data,dirs, sprintf('ROI %d: Raw data',sub));
                        end;
                        if ~isempty(ext(1,sub).adj)
                                % in this case, we only have one ext
                                % struct per ROI
                                extract_fig(ext(1,sub).adj,dirs, sprintf('ROI %d: Fitted data',sub));
                        end;
                end;
        end;
end;
spm('pointer','arrow');

function extract_fig(data,dirres,titlestr)

figure;
axp = axes;
hold on;
title(titlestr);

if ischar(dirres.dirsel)
  if strcmp(dirres.dirsel(1:3),'vox')
    data = data';
  end,
  switch dirres.dirsel(4:end)
   case 'no'
    l = plot(data);
   case 'avg'
    for k=1:size(data,1)
      plot(k,data(k,:),'rd');
    end;
    lmean=plot(mean(data,2));
  end;
  set(axp,'Xgrid','on', 'Xtick',1:size(data,1));
  ax = axis;
  axis([.5 size(data,1)+.5 ax(3:4)])
else
  % DTI plot
  for k=1:size(dirres.dirsel,1)
    cind=find(dirres.dirsel(k,:));
    cdata=data(cind,:);
    cmean(k)=mean(cdata(:));
    plot(k,cdata(:),'rd');
  end;

  lmean=plot(cmean);
  set(axp,'Xgrid','on', 'Xtick',1:size(cmean,2));
  ax = axis;
  axis([.5 size(cmean,2)+.5 ax(3:4)])
  
  if isfield(dirres,'b')
    lmeanb0=plot([1 size(dirres.dirsel,1)],cmean([1 1]),'r');
    set(axp,'Position',[.13 .31 .775 .615])
    
    axb=axes('Position',[.13 .21 .775 .05],...
	     'Xtick',.5:size(dirres.dirsel,1)+.5, 'XTickLabel',[],...
	     'Ytick',[], 'box','on');
    hold on;
    imagesc(dirres.b');
    colormap gray;
    axis tight;
    
    axd=axes('Position',[.13 .16 .775 .05],...
	     'Xtick',.5:size(dirres.dirsel,1)+.5, 'XTickLabel',[],...
	     'Ytick',[],'box','on');
    hold on;
    imagesc(dirres.g',[-1 1]);
    axis tight;
  end;
end;