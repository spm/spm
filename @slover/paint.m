function obj = paint(obj, params)
% Method to display slice overlay
% FORMAT obj = paint(obj, params)
%
% Inputs
% obj         - slice overlay object
% params      - optional structure containing extra display parameters
%               - refreshf - overrides refreshf in object
%               - clf      - overrides clf in object
%               - userdata - if 0, does not add object to userdata field
%               (see below)
%
% Outputs
% obj         - which may have been filled with defaults
%
% paint attaches the object used for painting to the 'UserData' field of
% the figure handle, unless instructed not to with 0 in userdata flag
%__________________________________________________________________________

% Matthew Brett
% $Id: paint.m 6623 2015-12-03 18:38:08Z guillaume $

fig_struct_fields = {'Position', 'Units'};

if nargin < 2
    params = [];
end
params = mars_struct('fillafromb', params, ...
    struct('refreshf', obj.refreshf,...
    'clf', obj.clf, ...
    'userdata', obj.userdata));

% Fill any needed defaults
obj = fill_defaults(obj);

% check if object can be painted
if isempty(obj.img)
    warning('slover:noImages', 'No images, object cannot be painted')
    return
end

% get coordinates for plane
X=1;Y=2;Z=3;
dims = obj.slicedef;
xmm = dims(X,1):dims(X,2):dims(X,3);
ymm = dims(Y,1):dims(Y,2):dims(Y,3);
zmm = obj.slices;
[y,x] = meshgrid(ymm,xmm');
vdims = [length(xmm),length(ymm),length(zmm)];

% no of slices, and panels (an extra for colorbars)
nslices = vdims(Z);
minnpanels = nslices;
cbars = 0;
if ~isempty(obj.cbar)
    cbars = length(obj.cbar);
    minnpanels = minnpanels+cbars;
end

% Get figure data.  The figure may be dead, in which case we may want to
% revive it.  If so, we set position etc as stored.
% If written to, the axes may be specified already
dead_f = ~ishandle(obj.figure);
figno = figure(obj.figure);
if dead_f
    set(figno, obj.figure_struct);
end

% (re)initialize axes and stuff

% check if the figure is set up correctly
if ~params.refreshf
    axisd = flipud(findobj(obj.figure, 'Type','axes','Tag', 'slice overlay panel'));
    npanels = length(axisd);
    if npanels < vdims(Z)+cbars;
        params.refreshf = 1;
    end
end
if params.refreshf
    % clear figure, axis store
    if params.clf, spm_clf(figno); end
    
    % prevent print inversion problems
    set(figno,'InvertHardCopy','off');
    
    % put copy of object into UserData for callbacks
    if params.userdata
        set(figno, 'UserData', obj);
    end
    
    % calculate area of display in pixels
    parea = obj.area.position;
    if ~strcmp(obj.area.units, 'pixels')
        ubu = get(obj.figure, 'units');
        set(obj.figure, 'units','pixels');
        tmp = get(obj.figure, 'Position');
        ascf = tmp(3:4);
        if ~strcmp(obj.area.units, 'normalized')
            set(obj.figure, 'units',obj.area.units);
            tmp = get(obj.figure, 'Position');
            ascf = ascf ./ tmp(3:4);
        end
        set(figno, 'Units', ubu);
        parea = parea .* repmat(ascf, 1, 2);
    end
    asz = parea(3:4);
    
    % by default, make most parsimonious fit to figure
    yxratio = length(ymm)*dims(Y,2)/(length(xmm)*dims(X,2));
    if isempty(obj.xslices)
        % iteration needed to optimize, surprisingly.  Thanks to Ian NS
        axlen(X,:)=asz(1):-1:1;
        axlen(Y,:)=yxratio*axlen(X,:);
        panels = floor(asz'*ones(1,size(axlen,2))./axlen);
        estnpanels = prod(panels);
        tmp = find(estnpanels >= minnpanels);
        if isempty(tmp)
            error('Whoops, cannot fit panels onto figure');
        end
        b = tmp(1); % best fitting scaling
        panels = panels(:,b);
        axlen = axlen(:, b);
    else
        % if xslices is specified, assume X is flush with X figure dimensions
        panels(X:Y,1) = [obj.xslices; 0];
        axlen(X:Y,1) = [asz(X)/panels(X); 0];
    end
    
    % Axis dimensions are in pixels.  This prevents aspect ratio rescaling
    panels(Y) = ceil(minnpanels/panels(X));
    axlen(Y) = axlen(X)*yxratio;
    
    % centre (etc) panels in display area as required
    divs = [Inf 2 1];the_ds = [0;0];
    the_ds(X) = divs(strcmp(obj.area.halign, {'left','center','right'}));
    the_ds(Y) = divs(strcmp(obj.area.valign, {'bottom','middle','top'}));
    startc = parea(1:2)' + (asz'-(axlen.*panels))./the_ds;
    
    % make axes for panels
    r=0;c=1;
    npanels = prod(panels);
    lastempty = npanels-cbars;
    axisd = nan(1, npanels);
    for i = 1:npanels
        % panel userdata
        if i<=nslices
            u.type = 'slice';
            u.no   = zmm(i);
        elseif i > lastempty
            u.type = 'cbar';
            u.no   = i - lastempty;
        else
            u.type = 'empty';
            u.no   = i - nslices;
        end
        axpos = [r*axlen(X)+startc(X) (panels(Y)-c)*axlen(Y)+startc(Y) axlen'];
        axisd(i) = axes(...
            'Parent',figno,...
            'XTick',[],...
            'XTickLabel',[],...
            'YTick',[],...
            'YTickLabel',[],...
            'Box','on',...
            'XLim',[1 vdims(X)],...
            'YLim',[1 vdims(Y)],...
            'Units', 'pixels',...
            'Position',axpos,...
            'Tag','slice overlay panel',...
            'UserData',u);
        r = r+1;
        if r >= panels(X)
            r = 0;
            c = c+1;
        end
    end
end

% sort out labels
if ischar(obj.labels)
    do_labels = ~strcmpi(obj.labels, 'none');
else
    do_labels = 1;
end
if do_labels
    labels = obj.labels;
    if iscell(labels.format)
        if length(labels.format)~=vdims(Z)
            error('Oh dear, expecting %d labels, but found %d', ...
                vdims(Z), length(labels.contents));
        end
    else
        % format string for mm from AC labelling
        fstr = labels.format;
        labels.format = cell(vdims(Z),1);
        acpt = obj.transform * [0 0 0 1]';
        for i = 1:vdims(Z)
            labels.format(i) = {sprintf(fstr,zmm(i)-acpt(Z))};
        end
    end
end

% split images into picture and contour
itypes = {obj.img(:).type};
tmp = strcmpi(itypes, 'contour');
contimgs = find(tmp);
pictimgs = find(~tmp);

% modify picture image colormaps with any new colours
npimgs = length(pictimgs);
lrn = zeros(npimgs,3);
cmaps = cell(npimgs);
for i = 1:npimgs
    cmaps(i)={obj.img(pictimgs(i)).cmap};
    lrnv = [obj.img(pictimgs(i)).outofrange, obj.img(pictimgs(i)).nancol];
    for j = 1:length(lrnv)
        if numel(lrnv{j})==1
            lrn(i,j) = lrnv{j};
        else
            cmaps(i) = {[cmaps{i}; lrnv{j}(1:3)]};
            lrn(i,j) = size(cmaps{i},1);
        end
    end
end

% cycle through slices displaying images
nvox = prod(vdims(1:2));
pandims = [vdims([2 1]) 3]; % NB XY transpose for display

zimg = zeros(pandims);
for i = 1:nslices
    ixyzmm = [x(:)';y(:)';ones(1,nvox)*zmm(i);ones(1,nvox)];
    img = zimg;
    for j = 1:npimgs
        thisimg = obj.img(pictimgs(j));
        i1 = sf_slice2panel(thisimg, ixyzmm, obj.transform, vdims);
        % rescale to colormap
        [csdata,badvals]= pr_scaletocmap(...
            i1,...
            thisimg.range(1),...
            thisimg.range(2),...
            thisimg.cmap,...
            lrn(j,:));
        % take indices from colormap to make true colour image
        iimg = reshape(cmaps{j}(csdata(:),:),pandims);
        tmp = repmat(logical(~badvals),[1 1 3]);
        if strcmpi(thisimg.type, 'truecolour')
            img(tmp) = img(tmp) + iimg(tmp)*thisimg.prop;
        else % split colormap effect
            img(tmp) = iimg(tmp)*thisimg.prop;
        end
    end
    % threshold out of range values
    img(img>1) = 1;
    
    image('Parent', axisd(i),...
        'ButtonDownFcn', obj.callback,...
        'CData',img);
    
    
    % do contour plot
    for j=1:length(contimgs)
        thisimg = obj.img(contimgs(j));
        i1 = sf_slice2panel(thisimg, ixyzmm, obj.transform, vdims);
        if any(any(isfinite(i1)))
            i1(i1<min(thisimg.range))=min(thisimg.range);
            i1(i1>max(thisimg.range))=max(thisimg.range);
            if ~any(diff(i1(isfinite(i1)))), continue, end % skip empty planes
            if mars_struct('isthere', thisimg, 'linespec')
                linespec = thisimg.linespec;
            else
                linespec = 'w-';
            end
            set(axisd(i),'NextPlot','add');
            if mars_struct('isthere', thisimg, 'contours')
                [c,h] = contour(axisd(i), i1, thisimg.contours, linespec);
            else
                [c,h] = contour(axisd(i), i1, linespec);
            end
            if ~isempty(h)
                if ~mars_struct('isthere', thisimg, 'linespec')
                    % need to reset colours; contour assumes you just set the figure's
                    % colormap, but to overlay coloured contours on a greyscale image,
                    % one needs to have everything in truecolour, setting individual
                    % contour lines to their appropriate colour.
                    % (updated for MATLAB 7 and above)
                    convals = get(h, 'LevelList');
                    if ~isempty(convals)
                        csdata = pr_scaletocmap(...
                            convals,...
                            thisimg.range(1),...
                            thisimg.range(2),...
                            thisimg.cmap,...
                            [1 size(thisimg.cmap,1) 1]);
                        colvals = thisimg.cmap(csdata(:),:)*thisimg.prop;
                        for ch = get(h, 'Children')' % (NB: transpose needed to loop)
                            CData = get(ch, 'CData'); % (CData is constant with final NaN)
                            colval = colvals(find(convals == CData(1), 1), :);
                            set(ch, 'EdgeColor', colval)
                        end
                    end
                end
                if mars_struct('isthere', thisimg, 'linewidth')
                    set(h, 'LineWidth', thisimg.linewidth);
                end
            end
        end
    end
    
    if do_labels
        text('Parent',axisd(i),...
            'Color', labels.colour,...
            'FontUnits', 'normalized',...
            'VerticalAlignment','bottom',...
            'HorizontalAlignment','left',...
            'Position', [1 1],...
            'FontSize',labels.size,...
            'ButtonDownFcn', obj.callback,...
            'String', labels.format{i});
    end
end
for i = (nslices+1):npanels
    set(axisd(i),'Color',[0 0 0]);
end
% add colorbar(s)
for i = 1:cbars
    axno = axisd(end-cbars+i);
    cbari = obj.img(obj.cbar(i));
    cml = size(cbari.cmap,1);
    p = get(axno, 'Position'); % position of last axis
    cw = p(3)*0.2;
    ch = p(4)*0.75;
    pc = p(3:4)/2;
    [axlims,idxs] = sort(cbari.range);
    a=axes(...
        'Parent',figno,...
        'XTick',[],...
        'XTickLabel',[],...
        'Units', 'pixels',...
        'YLim', axlims,...
        'FontUnits', 'normalized',...
        'FontSize', 0.075,...
        'YColor',[1 1 1],...
        'Tag', 'cbar',...
        'Box', 'off',...
        'Position',[p(1)+pc(1)-cw/2,p(2)+pc(2)-ch/2,cw,ch]...
        );
    image('Parent', a,...
        'YData', axlims(idxs),...
        'CData', reshape(cbari.cmap,[cml,1,3]));
end % colourbars

% Get stuff for figure, in case it dies later
obj.figure_struct = mars_struct('split', get(figno), fig_struct_fields);


%==========================================================================
function i1 = sf_slice2panel(img, xyzmm, transform, vdims)
% to voxel space of image
vixyz = (transform*img.vol.mat) \ xyzmm;
% raw data
if mars_struct('isthere', img.vol, 'imgdata')
    V = img.vol.imgdata;
else
    V = img.vol;
end
i1 = spm_sample_vol(V,vixyz(1,:),vixyz(2,:),vixyz(3,:), ...
    [img.hold img.background]);
if mars_struct('isthere', img, 'func')
    eval(img.func);
end
% transpose to reverse X and Y for figure
i1 = reshape(i1, vdims(1:2))';
