function [p,msg] = spm_eeval(str,Type,n,m)
% Expression evaluation
% FORMAT [p,msg] = spm_eeval(str,Type,n,m)
% Str  - Expression to work with
%
% Type - type of evaluation
%      - 's'tring
%      - 'e'valuated string
%          - 'n'atural numbers
%          - 'w'hole numbers
%          - 'i'ntegers
%          - 'r'eals
%      - 'c'ondition indicator vector
%
% n ('e', 'c' & 'p' types)
%          - Size of matrix requred
%          - NaN for 'e' type implies no checking - returns input as evaluated
%          - length of n(:) specifies dimension - elements specify size
%          - Inf implies no restriction
%          - Scalar n expanded to [n,1] (i.e. a column vector)
%            (except 'x' contrast type when it's [n,np] for np
%          - E.g: [n,1] & [1,n] (scalar n) prompt for an n-vector,
%                         returned as column or row vector respectively
%                 [1,Inf] & [Inf,1] prompt for a single vector,
%                         returned as column or row vector respectively
%                 [n,Inf] & [Inf,n] prompts for any number of n-vectors,
%                         returned with row/column dimension n respectively.
%                 [a,b] prompts for an 2D matrix with row dimension a and
%                         column dimension b
%                 [a,Inf,b] prompt for a 3D matrix with row dimension a,
%                         page dimension b, and any column dimension.
%          - 'c' type can only deal with single vectors
%          - NaN for 'c' type treated as Inf
%          - Defaults (missing or empty) to NaN
%
% m ('n', 'w', 'n1', 'w1', 'bn1' & 'bw1' types)
%          - Maximum value (inclusive)
%
% m ('r' type)
%          - Maximum and minimum values (inclusive)
%
% m ('c' type)
%       - Number of unique conditions required by 'c' type
%          - Inf implies no restriction
%          - Defaults (missing or empty) to Inf - no restriction
%
% p     - Result
%
% msg   - Explanation of why it didn't work
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_eeval.m 3756 2010-03-05 18:43:37Z guillaume $


if nargin<4, m=[]; end
if nargin<3, n=[]; end
if nargin<2, Type='e'; end
if nargin<1, str=''; end
if isempty(str), p='!'; msg='empty input'; return, end
switch lower(Type)
case 's'
    p = str; msg = '';
case 'e'
    p = evalin('base',['[',str,']'],'''!''');
    if ischar(p)
        msg = 'evaluation error';
    else
        [p,msg] = sf_SzChk(p,n);
    end
case 'n'
    p = evalin('base',['[',str,']'],'''!''');
    if ischar(p)
        msg = 'evaluation error';
    elseif any(floor(p(:))~=p(:)|p(:)<1) || ~isreal(p)
        p='!'; msg='natural number(s) required';
    elseif ~isempty(m) && any(p(:)>m)
        p='!'; msg=['max value is ',num2str(m)];
    else
        [p,msg] = sf_SzChk(p,n);
    end
case 'w'
    p = evalin('base',['[',str,']'],'''!''');
    if ischar(p)
        msg = 'evaluation error';
    elseif any(floor(p(:))~=p(:)|p(:)<0) || ~isreal(p)
        p='!'; msg='whole number(s) required';
    elseif ~isempty(m) && any(p(:)>m)
        p='!'; msg=['max value is ',num2str(m)];
    else
        [p,msg] = sf_SzChk(p,n);
    end
case 'i'
    p = evalin('base',['[',str,']'],'''!''');
    if ischar(p)
        msg = 'evaluation error';
    elseif any(floor(p(:))~=p(:)) || ~isreal(p)
        p='!'; msg='integer(s) required';
    else
        [p,msg] = sf_SzChk(p,n);
    end
case 'p'
    p = evalin('base',['[',str,']'],'''!''');
    if ischar(p)
        msg = 'evaluation error';
    elseif length(setxor(p(:)',m))
        p='!'; msg='invalid permutation';
    else
        [p,msg] = sf_SzChk(p,n);
    end
case 'r'
    p = evalin('base',['[',str,']'],'''!''');
    if ischar(p)
        msg = 'evaluation error';
    elseif ~isreal(p)
        p='!'; msg='real number(s) required';
    elseif ~isempty(m) && ( max(p)>max(m) || min(p)<min(m) )
        p='!'; msg=sprintf('real(s) in [%g,%g] required',min(m),max(m));
    else
        [p,msg] = sf_SzChk(p,n);
    end
case 'c'
    if isempty(m), m=Inf; end
    [p,msg] = icond(str,n,m);

otherwise
    error('unrecognised type');
end
return;

function [i,msg] = icond(i,n,m)
if nargin<3, m=Inf; end
if nargin<2, n=NaN; end
if any(isnan(n(:)))
    n=Inf;
elseif (length(n(:))==2 & ~any(n==1)) | length(n(:))>2
    error('condition input can only do vectors')
end
if nargin<2, i=''; end
if isempty(i), varargout={[],'empty input'}; return, end
msg = ''; i=i(:)';

if ischar(i)
    if i(1)=='0' & all(ismember(unique(i(:)),char(abs('0'):abs('9'))))
        %-Leading zeros in a digit list
        msg = sprintf('%s expanded',i);
        z = min(find([diff(i=='0'),1]));
        i = [zeros(1,z), icond(i(z+1:end))'];
    else
        %-Try an eval, for functions & string #s
        i = evalin('base',['[',i,']'],'i');
    end
end

if ischar(i)
    %-Evaluation error from above: see if it's an 'abab' or 'a b a b' type:
    [c,null,i] = unique(lower(i(~isspace(i))));
    if all(ismember(c,char(abs('a'):abs('z'))))
        %-Map characters a-z to 1-26, but let 'r' be zero (rest)
        tmp = c-'a'+1; tmp(tmp=='r'-'a'+1)=0;
        i   = tmp(i);
        msg = [sprintf('[%s] mapped to [',c),...
                sprintf('%d,',tmp(1:end-1)),...
                sprintf('%d',tmp(end)),']'];
    else
        i = '!'; msg = 'evaluation error';
    end
elseif ~all(floor(i(:))==i(:))
    i = '!'; msg = 'must be integers';
elseif length(i)==1 & prod(n)>1
    msg = sprintf('%d expanded',i);
    i = floor(i./10.^[floor(log10(i)+eps):-1:0]);
    i = i-[0,10*i(1:end-1)];
end

%-Check size of i & #conditions
if ~ischar(i), [i,msg] = sf_SzChk(i,n,msg); end
if ~ischar(i) & isfinite(m) & length(unique(i))~=m
    i = '!'; msg = sprintf('%d conditions required',m);
end
return;

function str = sf_SzStr(n,unused)
%=======================================================================
%-Size info string constuction
if nargin<2, l=0; else l=1; end
if nargin<1, error('insufficient arguments'); end;
if isempty(n), n=NaN; end
n=n(:); if length(n)==1, n=[n,1]; end; dn=length(n);
if any(isnan(n)) || (prod(n)==1 && dn<=2) || (dn==2 && min(n)==1 && isinf(max(n)))
    str = ''; lstr = '';
elseif dn==2 && min(n)==1
    str = sprintf('[%d]',max(n));   lstr = [str,'-vector'];
elseif dn==2 && sum(isinf(n))==1
    str = sprintf('[%d]',min(n));   lstr = [str,'-vector(s)'];
else
    str=''; for i = 1:dn
        if isfinite(n(i)), str = sprintf('%s,%d',str,n(i));
        else str = sprintf('%s,*',str); end
    end
    str = ['[',str(2:end),']'];     lstr = [str,'-matrix'];
end
if l, str=sprintf('\t%s',lstr); else str=[str,' ']; end


function [p,msg] = sf_SzChk(p,n,msg)
%=======================================================================
%-Size checking
if nargin<3, msg=''; end
if nargin<2, n=[]; end; if isempty(n), n=NaN; else n=n(:)'; end
if nargin<1, error('insufficient arguments'), end

if ischar(p) || any(isnan(n(:))), return, end
if length(n)==1, n=[n,1]; end

dn = length(n);
sp = size(p);

if dn==2 && min(n)==1
    %-[1,1], [1,n], [n,1], [1,Inf], [Inf,1] - vector - allow transpose
    %---------------------------------------------------------------
    i = min(find(n==max(n)));
    if n(i)==1 && max(sp)>1
        p='!'; msg='scalar required';
    elseif ndims(p)~=2 || ~any(sp==1) || ( isfinite(n(i)) && max(sp)~=n(i) )
        %-error: Not2D | not vector | not right length
        if isfinite(n(i)), str=sprintf('%d-',n(i)); else str=''; end
        p='!'; msg=[str,'vector required'];
    elseif sp(i)==1 && n(i)~=1
        p=p'; msg=[msg,' (input transposed)'];
    end

elseif dn==2 && sum(isinf(n))==1
    %-[n,Inf], [Inf,n] - n vector(s) required - allow transposing
    %---------------------------------------------------------------
    i = find(isfinite(n));
    if ndims(p)~=2 || ~any(sp==n(i))
        p='!'; msg=sprintf('%d-vector(s) required',min(n));
    elseif sp(i)~=n
        p=p'; msg=[msg,' (input transposed)'];
    end

else
    %-multi-dimensional matrix required - check dimensions
    %---------------------------------------------------------------
    if ndims(p)~=dn || ~all( size(p)==n | isinf(n) )
        p = '!'; msg='';
        for i = 1:dn
            if isfinite(n(i)), msg = sprintf('%s,%d',msg,n(i));
            else msg = sprintf('%s,*',msg); end
        end
        msg = ['[',msg(2:end),']-matrix required'];
    end
end
