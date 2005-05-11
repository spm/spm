function c = plus(a,b)

% FORMAT c = plus(a,b)
% 
% Simply overloads the 'plus' (+) function for datatype uint8.
% Beware that when values get over the range (255), it is set at 255.
%
% Written by c.phillips@ulg.ac.be, on 2004/11/25


sz_a = size(a); Nsz_a = length(sz_a);
sz_b = size(b); Nsz_b = length(sz_b);
if Nsz_a~=Nsz_b
    error('Matrices of different dimensions')
elseif ~all(sz_a==sz_b)
    error('Matrices of different dimensions')
else
    sz = sz_a; Nsz = Nsz_a;
end    

c = a;
switch Nsz
case 1
    tmp = double(a)+double(b);
    if tmp>255
        warning('Values went over the uint8 range!');
    end
    c = uint8(round(tmp));
case 2
    for ii=1:sz(2)
        tmp = double(a(:,ii))+double(b(:,ii));
        if any(tmp>255)
            warning('Values went over the uint8 range!');
        end
        c(:,ii) = uint8(round(tmp));
    end
case 3
    for ii=1:sz(3)
        tmp = double(a(:,:,ii))+double(b(:,:,ii));
        if any(tmp>255)
            warning('Values went over the uint8 range!');
        end
        c(:,:,ii) = uint8(round(tmp));
    end
otherwise
    error('this only works up to 3D matrices...');
end
