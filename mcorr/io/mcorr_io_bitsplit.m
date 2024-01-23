function d = mcorr_io_bitsplit(x,b,n)
% % bitsplit(X,B,N) splits the B-bit number X into signed N-bit array
% % X must be unsigned integer class
% % N ranges from 1 to B
% % B is a multiple of N
%
% sign = repmat((b:-n:n)',1,size(x,1));
% x = repmat(x',b/n,1);
% d = double(bitand(bitshift(x,flipud(sign-b)),2^n-1)) ...
%   - double(bitget(x,sign))*2^n;
