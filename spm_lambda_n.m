function l_n = spm_lambda_n(n)
% Computation of K.W.'s correction factor lambda_n. Integral of sjk.m 
% Looks up values in a precomputed table (actual computation time is some 
% minutes)
% 
% n : degrees of freedom + 1;
% l_n : see JCBF (K.J. Worsely, pp. 900 - 918, 1992).
%_______________________________________________________________________
% %W% Stefan Kiebel %E%

l =  [inf    inf    inf   1.7595 1.4351 1.3036 1.2333 1.1894 1.1593 1.1375,...
     1.1209 1.1079 1.0974 1.0888 1.0815 1.0754 1.0701 1.0655 1.0615 1.0579,...
     1.0548 1.0519 1.0494 1.0470 1.0449 1.0430 1.0412 1.0396 1.0381 1.0367,...
     1.0354 1.0342 1.0331 1.0320 1.0310 1.0301 1.0292 1.0284 1.0276 1.0268,...
     1.0261 1.0255 1.0248 1.0242 1.0237 1.0231 1.0226 1.0221 ];
l_n = l(n);

