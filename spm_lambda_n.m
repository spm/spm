function l_n = spm_lambda_n(n,quad)
% FORMAT l_n = spm_lambda_n(n)
%
% n    -   df + 1 : Degrees of freedon of the t-field, plus one
% quad -   Force use numerical integration rather than lookup table
%_______________________________________________________________________
%
% Returns the ratio of the Variance-Covariance matrix of partial
% derivatives of a Gaussianised t-field to the components of the
% t-field.
%
% \lambda_n is the integral of spm_lambda_n_int over all real values.
%
% The formula is that on p917 of Worsley et al. (1992) Worsley et al.,
% (1992) "A Three Dimensional Statistical Analysis for CBF Activation
% Studies in Human Brain", Journal of Cerebral Blood and Metabolism
% 12:900-918
%
% A precomputed lookup table is used for low n, since actual
% computation time can be tens of seconds. Working for obtaining the
% limits for the quadrature are presented as comments at the end of the 
% code.
%
%_______________________________________________________________________
% %W% Stefan Kiebel, Andrew Holmes %E%

% Check parameters
%-----------------------------------------------------------------------
if (nargin==0), error('Specify n'), end
if ( (n<1) | (floor(n)~=ceil(n)) )
	error('n must be a positive integer'), end
if (nargin<2), quad=0; end


% Define lookup table & quadrature parameters
%=======================================================================
lookup =  [
inf    inf    inf    1.7580 1.4339 1.3036 1.2333 1.1894 1.1593 1.1375,...
1.1209 1.1079 1.0974 1.0888 1.0815 1.0754 1.0701 1.0655 1.0615 1.0579,...
1.0548 1.0519 1.0494 1.0470 1.0449 1.0430 1.0412 1.0396 1.0381 1.0367,...
1.0354 1.0342 1.0331 1.0320 1.0310 1.0301 1.0292 1.0284 1.0276 1.0268,...
1.0261 1.0255 1.0248 1.0242 1.0237 1.0231 1.0226 1.0221 ];

tol = 1e-6;
nLim = [   4,   5,   6, 10, 15, 25, 30, 50,  100, inf];
tLim = [4000, 500, 100, 20, 15, 10, 10,  7.5,  5,   5];


% Computation
%=======================================================================
if ( n <= length(lookup) ) & ~quad
	% Select value from lookup table (if possible)
	%---------------------------------------------------------------
	l_n = lookup(n);
else
	% Quadrature
	%---------------------------------------------------------------
	if any([1,2,3]==n), l_n = inf; return, end
	iLim = min(find(n<=nLim));
	Lim = tLim(iLim);
	l_n = 2*quad8('spm_lambda_n_int',0,Lim,tol,[],n);
end


return

%   R O U G H   W O R K
%=======================================================================
% % n=4
% %-----------------------------------------------------------------------
% 2 *quad8('spm_lambda_n_int',0,15,1e-6,1,4) % = 1.72142839782810
% quad8('spm_lambda_n_int',-15,15,1e-6,1,4) % = 1.72142839782810
% 
% a = quad8('spm_lambda_n_int',0,10,1e-6,1,4) % = 0.84763748822928
% b = quad8('spm_lambda_n_int',10,50,1e-6,1,4) % = 0.02741516955612
% c = quad8('spm_lambda_n_int',50,5000,1e-6,1,4) % = 0.00392945357381
% (a+b+c)*2 % = 1.75796422271842
% c = quad8('spm_lambda_n_int',50,4000,1e-6,1,4) % = 0.00392411503994
% (a+b+c)*2 % = 1.75795354565068
% 
% 2 * quad8('spm_lambda_n_int',5000,6000,1e-6,1,4) % = 6.938622917610468e-06
% 2 * quad8('spm_lambda_n_int',6000,8000,1e-6,1,4) % = 8.433424081192626e-06
% 2 * quad8('spm_lambda_n_int',8000,10000,1e-6,1,4) % = 4.907710729890348e-06
% 2 * quad8('spm_lambda_n_int',10000,12000,1e-6,1,4) % = 3.1957e-06
% %-Strange rounding/underflow effects?
% 
% 2 * quad8('spm_lambda_n_int',4000,5000,1e-6,1,4) % = 1.067704696790532e-05
% 2 * quad8('spm_lambda_n_int',3000,4000,1e-6,1,4) % = 1.839671355433203e-05
% 2 * quad8('spm_lambda_n_int',2000,4000,1e-6,1,4) % = 5.697627206413217e-05
% 2 * quad8('spm_lambda_n_int',1000,4000,1e-6,1,4) % = 1.826777256412703e-04
% 
% % l_n = 1.7580 to 4 d.p.
% 
% 2 *quad8('spm_lambda_n_int',0,5000,1e-6,1,4) % = 1.75793784191677
% 2 *quad8('spm_lambda_n_int',0,5500,1e-6,1,4) % = 1.75791892157342 % ? smaller!
% 2 *quad8('spm_lambda_n_int',0,4000,1e-6,1,4) % = 1.75795063335648 % OK
% 2 *quad8('spm_lambda_n_int',0,3000,1e-6,1,4) % = 1.75793531629395 % Nah
% 
% % tLim = 4000
% 
% % n=5
% %-----------------------------------------------------------------------
% a = quad8('spm_lambda_n_int',0,10,1e-6,1,5) % = 0.71241952430330
% b = quad8('spm_lambda_n_int',10,50,1e-6,1,5) % = 0.00442671986537
% c = quad8('spm_lambda_n_int',50,4000,1e-6,1,5) % = 1.087591995552818e-04
% (a+b+c)*2 % = 1.43391000673645
% 
% 2 * quad8('spm_lambda_n_int',4000,5000,1e-6,1,5) % = 5.936265150284404e-09
% 
% 2 * quad8('spm_lambda_n_int',3000,4000,1e-6,1,5) % = 1.320775109485138e-08
% 2 * quad8('spm_lambda_n_int',2000,3000,1e-6,1,5) % = 3.962231338604665e-08
% 2 * quad8('spm_lambda_n_int',1000,2000,1e-6,1,5) % = 2.337274478980690e-07
% 2 * quad8('spm_lambda_n_int',500,2000,1e-6,1,5) % = 1.279518377610225e-06
% 2 * quad8('spm_lambda_n_int',250,2000,1e-6,1,5) % = 6.024028702195433e-06
% 
% % = 1.4339 to 4 d.p.
% 
% 2 *quad8('spm_lambda_n_int',0,2000,1e-6,1,5) % = 1.43390995861575 %
% 2 *quad8('spm_lambda_n_int',0,1000,1e-6,1,5) % = 1.43390972023890 %
% 2 *quad8('spm_lambda_n_int',0,500,1e-6,1,5) % = 1.43390867444656 % v/
% 2 *quad8('spm_lambda_n_int',0,250,1e-6,1,5) % = 1.43390392993730 %
% 
% % tLim = 500
% 
% % n=6
% %-----------------------------------------------------------------------
% a = quad8('spm_lambda_n_int',0,10,1e-6,1,6) % = 0.65083961564008
% b = quad8('spm_lambda_n_int',10,50,1e-6,1,6) % = 9.696146070187920e-04
% c = quad8('spm_lambda_n_int',50,500,1e-6,1,6) % = 4.601187004685710e-06
% (a+b+c)*2 % = 1.30362766286820
% 
% 2 * quad8('spm_lambda_n_int',500,1000,1e-6,1,6) % = 4.941646644133863e-09
% 
% 2 * quad8('spm_lambda_n_int',250,500,1e-6,1,6) % = 4.499785662841468e-08
% 2 * quad8('spm_lambda_n_int',200,500,1e-6,1,6) % = 9.771267409599348e-08
% 2 * quad8('spm_lambda_n_int',100,500,1e-6,1,6) % = 9.572244191990791e-07
% 2 * quad8('spm_lambda_n_int',50,500,1e-6,1,6) % = 9.202374009371421e-06
% 
% % = 1.3036 to 4 d.p.
% 
% 2 *quad8('spm_lambda_n_int',0,100,1e-6,1,6) % = 1.30362670564138
% 2 *quad8('spm_lambda_n_int',0,250,1e-6,1,6) % = 
% 
% % tLim = 100
% 
% % n=7
% %-----------------------------------------------------------------------
% a = quad8('spm_lambda_n_int',0,10,1e-6,1,7) % = 0.61640412951836
% b = quad8('spm_lambda_n_int',10,50,1e-6,1,7) % = 2.597643087526469e-04
% c = quad8('spm_lambda_n_int',50,500,1e-6,1,7) % = 2.469884368055823e-07
% (a+b+c)*2 % = 1.23332828163110
% 
% 2 * quad8('spm_lambda_n_int',500,600,1e-6,1,7) % = 1.572271187812036e-11
% 
% 2 * quad8('spm_lambda_n_int',250,500,1e-6,1,7) % = 5.060745490709128e-10
% 2 * quad8('spm_lambda_n_int',100,500,1e-6,1,7) % = 2.562067732151013e-08
% 2 * quad8('spm_lambda_n_int',50,500,1e-6,1,7) % = 4.939768736111646e-07
% 2 * quad8('spm_lambda_n_int',25,500,1e-6,1,7) % = 9.814458415082894e-06
% 2 * quad8('spm_lambda_n_int',15,500,1e-6,1,7) % = 9.012865417887753e-05
% 
% % = 1.2333 to 4 d.p.
% 
% 2 *quad8('spm_lambda_n_int',0,50,1e-6,1,7) % = 1.23332778765327
% 2 *quad8('spm_lambda_n_int',0,25,1e-6,1,7) % = 1.23331846718079
% 2 *quad8('spm_lambda_n_int',0,20,1e-6,1,7) % = 1.23330244328157
% 
% % tLim = 20
% 
% % n=8
% %-----------------------------------------------------------------------
% a = quad8('spm_lambda_n_int',0,5,1e-6,1,8) % = 0.59217316595614
% b = quad8('spm_lambda_n_int',5,10,1e-6,1,8) % = 0.00243802709776
% c = quad8('spm_lambda_n_int',10,100,1e-6,1,8) % = 8.033062363764178e-05
% (a+b+c)*2 % = 1.18938304735507
% 
% 2 * quad8('spm_lambda_n_int',100,200,1e-6,1,8) % = 7.832894955025816e-10
% 2 * quad8('spm_lambda_n_int',100,300,1e-6,1,8) % = 8.021950290943954e-10
% 
% 2 * quad8('spm_lambda_n_int',25,100,1e-6,1,8) % = 1.240917219635505e-06
% 2 * quad8('spm_lambda_n_int',20,100,1e-6,1,8) % = 4.082219546853740e-06
% 2 * quad8('spm_lambda_n_int',15,100,1e-6,1,8) % = 1.890715810307602e-05
% 2 * quad8('spm_lambda_n_int',10,100,1e-6,1,8) % = 1.606612472752836e-04
% 
% % = 1.1894 to 4 d.p.
% 
% 2 *quad8('spm_lambda_n_int',0,25,1e-6,1,8) % = 1.18938180644381
% 2 *quad8('spm_lambda_n_int',0,20,1e-6,1,8) % = 1.18937896513552
% 2 *quad8('spm_lambda_n_int',0,15,1e-6,1,8) % = 1.18936414023076
% 
% % tLim = 20
% 
% 
% % n=11
% %-----------------------------------------------------------------------
% a = quad8('spm_lambda_n_int',0,5,1e-6,1,11) % = 0.55979886803584
% b = quad8('spm_lambda_n_int',5,10,1e-6,1,11) % = 6.521051049054471e-04
% c = quad8('spm_lambda_n_int',10,50,1e-6,1,11) % = 4.188008482204222e-06
% (a+b+c)*2 % = 1.12091032229846
% 
% 2 * quad8('spm_lambda_n_int',50,100,1e-6,1,11) % = 1.483314099751780e-11
% 
% 2 * quad8('spm_lambda_n_int',25,50,1e-6,1,11) % = 4.702617176780218e-09
% 2 * quad8('spm_lambda_n_int',20,50,1e-6,1,11) % = 2.996826279052569e-08
% 2 * quad8('spm_lambda_n_int',15,50,1e-6,1,11) % = 3.200078929834791e-07
% 2 * quad8('spm_lambda_n_int',10,50,1e-6,1,11) % = 8.376016964408444e-06
% 
% % = 1.1209 to 4 d.p.
% 
% 2 *quad8('spm_lambda_n_int',0,25,1e-6,1,11) % = 1.12091031759683
% 2 *quad8('spm_lambda_n_int',0,20,1e-6,1,11) % = 1.12091029232985
% 2 *quad8('spm_lambda_n_int',0,15,1e-6,1,11) % = 1.12091000228971
% 
% % tLim = 15
% 
% % n=16
% %-----------------------------------------------------------------------
% a = quad8('spm_lambda_n_int',0,5,1e-6,1,16) % = 0.53754543155844
% b = quad8('spm_lambda_n_int',5,10,1e-6,1,16) % = 1.523493195225369e-04
% c = quad8('spm_lambda_n_int',10,25,1e-6,1,16) % = 9.455380698933340e-08
% (a+b+c)*2 % = 1.07539575086354
% 
% 2 * quad8('spm_lambda_n_int',25,35,1e-6,1,16) % = 1.583498755148690e-12
% 
% 2 * quad8('spm_lambda_n_int',20,25,1e-6,1,16) % = 2.824837382736746e-11
% 2 * quad8('spm_lambda_n_int',15,25,1e-6,1,16) % = 1.226664801743324e-09
% 2 * quad8('spm_lambda_n_int',10,25,1e-6,1,16) % = 1.891076139786668e-07
% 2 * quad8('spm_lambda_n_int',5,25,1e-6,1,16) % = 3.048877466600576e-04
% 
% % = 1.0754 to 4 d.p.
% 
% 2 *quad8('spm_lambda_n_int',0,25,1e-6,1,16) % = 1.07539574325938
% 2 *quad8('spm_lambda_n_int',0,20,1e-6,1,16) % = 1.07539574323426
% 2 *quad8('spm_lambda_n_int',0,15,1e-6,1,16) % = 1.07539574203510
% 2 *quad8('spm_lambda_n_int',0,10,1e-6,1,16) % = 1.07539574203510
% 
% % tLim = 10
% 
% % n=26
% %-----------------------------------------------------------------------
% a = quad8('spm_lambda_n_int',0,5,1e-6,1,26) % = 0.52146933383945
% b = quad8('spm_lambda_n_int',5,10,1e-6,1,26) % = 2.881927073548046e-05
% c = quad8('spm_lambda_n_int',10,15,1e-6,1,26) % = 4.317857591419682e-10
% (a+b+c)*2 % = 1.04299630708394
% 
% 2 * quad8('spm_lambda_n_int',15,16,1e-6,1,26) % = 1.682442368240151e-13
% 
% 2 * quad8('spm_lambda_n_int',10,15,1e-6,1,26) % = 8.635715182839365e-10
% 2 * quad8('spm_lambda_n_int',5,15,1e-6,1,26) % = 5.763940504250183e-05
% 
% % = 1.0430 to 4 d.p.
% 
% 2 *quad8('spm_lambda_n_int',0,15,1e-6,1,26) % = 1.04299630708394
% 2 *quad8('spm_lambda_n_int',0,10,1e-6,1,26) % = 1.04299630622037
% 
% % tLim = 10
% 
% % n=31
% %-----------------------------------------------------------------------
% a = quad8('spm_lambda_n_int',0,5,1e-6,1,31) % = 0.51767811649975
% b = quad8('spm_lambda_n_int',5,10,1e-6,1,31) % = 1.694769624816858e-05
% (a+b)*2 % = 1.03539012839200
% 
% 2 * quad8('spm_lambda_n_int',10,15,1e-6,1,31) % = 1.110384094791734e-10
% 
% 2 * quad8('spm_lambda_n_int',8,10,1e-6,1,31) % = 1.234623962986153e-08
% 2 * quad8('spm_lambda_n_int',7,10,1e-6,1,31) % = 1.588658961595441e-07
% 2 * quad8('spm_lambda_n_int',6,10,1e-6,1,31) % = 2.250230282304706e-06
% 2 * quad8('spm_lambda_n_int',5,10,1e-6,1,31) % = 3.389539249633716e-05
% 
% % = 1.0354 to 4 d.p.
% 
% 2 *quad8('spm_lambda_n_int',0,15,1e-6,1,31) % = 1.03539012850303
% 2 *quad8('spm_lambda_n_int',0,10,1e-6,1,31) % = 1.03539012839200
% 2 *quad8('spm_lambda_n_int',0,8,1e-6,1,31) % = 1.03539011604575
% 2 *quad8('spm_lambda_n_int',0,7.5,1e-6,1,31) % = 1.03539008464327
% 
% % tLim = 7.5
% 
% 
% % n=51
% %-----------------------------------------------------------------------
% a = quad8('spm_lambda_n_int',0,5,1e-6,1,51) % = 0.51035713830147
% b = quad8('spm_lambda_n_int',5,10,1e-6,1,51) % = 4.733784538394881e-06
% (a+b)*2 % = 1.02072374417201
% 
% 2 * quad8('spm_lambda_n_int',10,12,1e-6,1,51) % = 3.035677420930937e-13
% 
% 2 * quad8('spm_lambda_n_int',8,10,1e-6,1,51) % = 2.672575234920359e-10
% 2 * quad8('spm_lambda_n_int',7,10,1e-6,1,51) % = 8.908917552312104e-09
% 2 * quad8('spm_lambda_n_int',6,10,1e-6,1,51) % = 3.004541512640518e-07
% 2 * quad8('spm_lambda_n_int',5,10,1e-6,1,51) % = 9.467569076789761e-06
% 
% % = 1.0207 to 4 d.p.
% 
% 2 *quad8('spm_lambda_n_int',0,8,1e-6,1,51) % = 1.02072374390475
% 2 *quad8('spm_lambda_n_int',0,7,1e-6,1,51) % = 1.02072373526458
% 2 *quad8('spm_lambda_n_int',0,6,1e-6,1,51) % = 1.02072344371786
% 2 *quad8('spm_lambda_n_int',0,5,1e-6,1,51) % = 1.02071427660293
% 
% % tLim = 5
