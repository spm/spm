function c = times(a,b)

% FORMAT c = times(a,b)
% 
% Simply overloads the 'times' (*) function for datatype uint8.
% Only the multiplication of a matrix by a scalar is overloaded !
% Beware that when values get over the range (255), it is set at 255.
%
% Written by c.phillips@ulg.ac.be, on 2004/11/25


sz = size(a); Nsz = length(sz);
if prod(size(b))>1
    error('Only multiply a matrix by a scalar')
end    

c = a;
switch Nsz
case 1
    tmp = double(a)*double(b);
    if tmp>255
        warning('Values went over the uint8 range!');
    end
    c = uint8(round(tmp));
case 2
    for ii=1:sz(2)
        tmp = double(a(:,ii))*double(b);
        if any(tmp>255)
            warning('Values went over the uint8 range!');
        end
        c(:,ii) = uint8(round(tmp));
    end
case 3
    for ii=1:sz(3)
        tmp = double(a(:,:,ii))*double(b);
        if any(tmp>255)
            warning('Values went over the uint8 range!');
        end
        c(:,:,ii) = uint8(round(tmp));
    end
otherwise
    error('This only works with matrices up to 3D...');
end
