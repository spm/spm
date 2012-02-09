function y = identity(d)

y = zeros([d(1:3) 3]);
[y(:,:,:,1),y(:,:,:,2),y(:,:,:,3)] = ndgrid(single(1:d(1)),single(1:d(2)),single(1:d(3)));
