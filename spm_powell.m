function [p,f] = spm_powell(p,xi,tolsc,func,varargin)
% Powell optimisation method
% FORMAT [p,f] = spm_powell(p,xi,tolsc,func,varargin)
% 	p        - Starting parameter values
%	xi       - columns containing directions in which to begin
% 	           searching.
% 	tolsc    - stopping criteria
% 	           - optimisation stops when
% 	             sqrt(sum(((p-p_prev)./tolsc).^2))<1
% 	func     - name of evaluated function
% 	varargin - remaining arguments to func (after p)
%
% 	p        - final parameter estimates
% 	f        - function value at minimum
%
%_______________________________________________________________________
% Method is based on Powell's optimisation method from Numerical Recipes
% in C (Press, Flannery, Teukolsky & Vetterling).
%_______________________________________________________________________
% %W% John Ashburner %E%

p     = p(:);
f     = feval(func,p,varargin{:});
for iter=1:128,
	fprintf('iteration %d...\n', iter);
	pp   = p;
	fp   = f;
	del  = 0;
	for i=1:length(p),
		ft = f;
		[p,junk,f] = linmin(p,xi(:,i),func,f,tolsc,varargin{:});
		if abs(ft-f) > del,
			del  = abs(ft-f);
			ibig = i;
		end;
	end;
	if sqrt(sum(((p(:)-pp(:))./tolsc(:)).^2))<1, return; end;
	ft = feval(func,2.0*p-pp,varargin{:});
	if ft < f,
		[p,xi(:,ibig),f] = linmin(p,p-pp,func,f,tolsc,varargin{:});
	end;
end;
warning('Too many optimisation iterations');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [p,pi,f] = linmin(p,pi,func,f,tolsc,varargin)
% Line search for minimum.

global lnm % used in linmineval
lnm      = struct('p',p,'pi',pi,'func',func,'args',[]);
lnm.args = varargin;

linmin_plot('Init', 'Line Minimisation','Function','Parameter Value');
linmin_plot('Set', 0, f);

tol      = 1/sqrt(sum((pi(:)./tolsc(:)).^2));
t        = bracket(f);
[f,pmin] = brents(t,tol);
pi       = pi*pmin;
p        = p + pi;

for i=1:length(p), fprintf('%-8.4g ', p(i)); end;
fprintf('| %.5g\n', f);
linmin_plot('Clear');

return;
%_______________________________________________________________________

%_______________________________________________________________________
function f = linmineval(p)
% Reconstruct parameters and evaluate.

global lnm % defined in linmin
pt = lnm.p+p.*lnm.pi;
f  = feval(lnm.func,pt,lnm.args{:});
linmin_plot('Set',p,f);
return;
%_______________________________________________________________________

%_______________________________________________________________________
function t = bracket(f)
% Bracket the minimum (t(2)) between t(1) and t(2)

gold   = (1+sqrt(5))/2; % Golden ratio

t(1)   = struct('p',0,'f',f);
t(2).p = 1;
t(2).f = linmineval(t(2).p);

% if not better then swap
if t(2).f > t(1).f,
	tmp  = t(1);
	t(1) = t(2);
	t(2) = tmp;
end;

t(3).p = t(2).p + gold*(t(2).p-t(1).p);
t(3).f = linmineval(t(3).p);

while t(2).f > t(3).f,

	% fit a polynomial to t
	tmp = cat(1,t.p)-t(2).p;
	pol = pinv([ones(3,1) tmp tmp.^2])*cat(1,t.f);

	% minimum is when gradient of polynomial is zero
	% sign of pol(3) (the 2nd deriv) should be +ve
	d   = -pol(2)/(2*pol(3)+eps);
	u.p = t(2).p+d;

	ulim = t(2).p+100*(t(3).p-t(2).p);
	if (t(2).p-u.p)*(u.p-t(3).p) > 0.0,
		% u is between t(2) and t(3)
		u.f = linmineval(u.p);
		if u.f < t(3).f,
			% minimum between t(2) and t(3) - done
			t(1) = t(2);
			t(2) = u;
			return;
		elseif u.f > t(2).f,
			% minimum between t(1) and u - done
			t(3) = u;
			return;
		end;
		% try golden search instead
		u.p = t(3).p+gold*(t(3).p-t(2).p);
		u.f = linmineval(u.p);

	elseif (t(3).p-u.p)*(u.p-ulim) > 0.0
		% u is between t(3) and ulim
		u.f = linmineval(u.p);
		if u.f < t(3).f,
			% still no minimum as function is still decreasing
			% t(1) = t(2);
			t(2) = t(3);
			t(3) = u;
			u.p  = t(3).p+gold*(t(3).p-t(2).p);
			u.f  = linmineval(u.p);
		end;

	elseif (u.p-ulim)*(ulim-t(3).p) >= 0.0,
		% gone too far - constrain it
		u.p = ulim;
		u.f = linmineval(u.p);

	else,
		% try golden search instead
		u.p = t(3).p+gold*(t(3).p-t(2).p);
		u.f = linmineval(u.p);
	end;

	% Move all 3 points along
	t(1) = t(2);
	t(2) = t(3);
	t(3) = u;
end;
return;
%_______________________________________________________________________

%_______________________________________________________________________
function [f,p] = brents(t, tol)
% Brent's method for line searching - given that minimum is bracketed

% 1 - golden ratio
Cgold = 1 - (sqrt(5)-1)/2;

% Current and previous displacements
d     = Inf;
pd    = Inf;

% t(1) and t(3) bracket the minimum
if t(1).p>t(3).p,
	brk(1) = t(3).p;
	brk(2) = t(1).p;
else,
	brk(1) = t(1).p;
	brk(2) = t(3).p;
end;

% sort t into best first order
tmp  = t(1);
t(1) = t(2);
t(2) = tmp;
if t(2).f>t(3).f,
	tmp  = t(2);
	t(2) = t(3);
	t(3) = tmp;
end;

for iter=1:128,
	% check stopping criterion
	if abs(t(1).p - 0.5*(brk(1)+brk(2)))+0.5*(brk(2)-brk(1)) <= 2*tol,
		p = t(1).p;
		f = t(1).f;
		return;
	end;

	% keep last two displacents
	ppd = pd;
	pd  = d;

	% fit a polynomial to t
	tmp = cat(1,t.p)-t(1).p;
	pol = pinv([ones(3,1) tmp tmp.^2])*cat(1,t.f);

	% minimum is when gradient of polynomial is zero
	d   = -pol(2)/(2*pol(3)+eps);
	u.p = t(1).p+d;

	% check so that displacement is less than the last but two,
	% that the displaced point is between the brackets
	% and (not sure if it is necessary) that the solution is a minimum
	% rather than a maximum
	eps2 = 2*eps*abs(t(1).p)+eps;
	if abs(d) >= abs(ppd)/2 | u.p <= brk(1)+eps2 | u.p >= brk(2)-eps2 | pol(3)<=0,
		% if criteria are not met, then golden search into the larger part
		if t(1).p >= 0.5*(brk(1)+brk(2)),
			d = Cgold*(brk(1)-t(1).p);
		else,
			d = Cgold*(brk(2)-t(1).p);
		end;
		u.p = t(1).p+d;
	end;

	% FUNCTION EVALUATION
	u.f = linmineval(u.p);

	% Insert the new point into the appropriate position and update
	% the brackets if necessary
	if u.f <= t(1).f,
		if u.p >= t(1).p, brk(1)=t(1).p; else, brk(2)=t(1).p; end;
		t(3) = t(2);
		t(2) = t(1);
		t(1) = u;
	else,
		if u.p < t(1).p, brk(1)=u.p; else, brk(2)=u.p; end;
		if u.f <= t(2).f | t(1).p==t(2).p,
			t(3) = t(2);
			t(2) = u;
		elseif u.f <= t(3).f | t(1).p==t(3).p | t(2).p==t(3).p,
			t(3) = u;
		end;
	end;
end;
warning('Too many iterations in Brents');
return;
%_______________________________________________________________________

%_______________________________________________________________________
function linmin_plot(action,arg1,arg2,arg3,arg4)
% Visual output for line minimisation
global linminplot
%-----------------------------------------------------------------------
if (nargin == 0)
	linmin_plot('Init');
else
	% initialize
	%---------------------------------------------------------------
	if (strcmp(lower(action),'init'))
		if (nargin<4)
			arg3 = 'Function';
			if (nargin<3)
				arg2 = 'Value';
				if (nargin<2)
					arg1 = 'Line minimisation';
				end
			end
		end
		fg = spm_figure('FindWin','Interactive');

		if ~isempty(fg)
			linminplot = struct('pointer',get(fg,'Pointer'),'name',get(fg,'Name'),'ax',[]);
			linmin_plot('Clear');
			set(fg,'Pointer','watch');
			% set(fg,'Name',arg1);
			linminplot.ax = axes('Position', [0.15 0.1 0.8 0.75],...
				'Box', 'on','Parent',fg);
			lab = get(linminplot.ax,'Xlabel');
			set(lab,'string',arg3,'FontSize',10);
			lab = get(linminplot.ax,'Ylabel');
			set(lab,'string',arg2,'FontSize',10);
			lab = get(linminplot.ax,'Title');
			set(lab,'string',arg1);
			line('Xdata',[], 'Ydata',[],...
				'LineWidth',2,'Tag','LinMinPlot','Parent',linminplot.ax,...
				'LineStyle','-','Marker','o');
			drawnow;
		end

	% reset
	%---------------------------------------------------------------
	elseif (strcmp(lower(action),'set'))
		F = spm_figure('FindWin','Interactive');
		br = findobj(F,'Tag','LinMinPlot');
		if (~isempty(br))
			[xd,indx] = sort([get(br,'Xdata') arg1]);
			yd = [get(br,'Ydata') arg2];
			yd = yd(indx);
			set(br,'Ydata',yd,'Xdata',xd);
			drawnow;
		end

	% clear
	%---------------------------------------------------------------
	elseif (strcmp(lower(action),'clear'))
		fg = spm_figure('FindWin','Interactive');
		if isstruct(linminplot),
			if ishandle(linminplot.ax), delete(linminplot.ax); end;
			set(fg,'Pointer',linminplot.pointer);
			set(fg,'Name',linminplot.name);
		end;
		spm_figure('Clear',fg);
		drawnow;
	end;
end
%_______________________________________________________________________
