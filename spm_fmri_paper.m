
% plots individual transients following AnCova by using variables in working
% memory extant after calling spm_plot
%
%__________________________________________________________________________
% %W% %E%

D     = [1 51 71 91;11 31 61 101;21 41 81 111] - 1;
%----------------------------------------------------------------------------
figure(3); spm_clf
hold on
x     = XA(:,i);
% x     = C*(pinv(C)*x);
x     = [x; x];
t     = [0:11]*RT;
for j = 1:size(D,2)
	d  = [1:12] + D(1,j);
	d  = x(d); % d = d - d(2);
	plot(t,d,':')
end
for j = 1:size(D,2)
	d  = [1:12] + D(2,j);
	d  = x(d); % d = d - d(2);
	plot(t,d,'-.')
end
for j = 1:size(D,2)
	d  = [1:12] + D(3,j);
	d  = x(d); % d = d - d(2);
	plot(t,d,'-')
end
hold off

%----------------------------------------------------------------------------
title(sprintf('Hemodynamic responses at %0.0i, %0.0i, %0.0i mm',L),'FontSize',16)
xlabel('time {seconds}')
ylabel('adjusted activity')
axis square
grid on

