function varargout = spm_axis(varargin)
% AXIS  Control axis scaling and appearance.

if nargout
    [varargout{1:nargout}] = axis(varargin{:});
else
    try
        axis(varargin{:});
    end
end

if nargin == 1 && any(strcmpi(varargin{1},{'tight','scale'}))
    spm_axis(gca,varargin{1});
elseif nargin == 2 && allAxes(varargin{1}) && strcmpi(varargin{2},'tight')
    for i = 1:numel(varargin{1})
        lm = get(varargin{1}(i),'ylim');
        set(varargin{1}(i),'ylim',lm + [-1 1]*diff(lm)/16);
    end
elseif nargin == 2 && allAxes(varargin{1}) && strcmpi(varargin{2},'scale')
    for i = 1:numel(varargin{1})
        lm = get(varargin{1}(i),'ylim');
        set(varargin{1}(i),'ylim',[0 lm(2)*(1 + 1/16)]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = allAxes(h)

result = all(ishghandle(h)) && ...
         length(findobj(h,'type','axes','-depth',0)) == length(h);
