function d = mcorr_io_bitsign(x,n)
% bitsign(X,N) returns signed double value from unsigned N-bit number X.
% This is equivalent to bitsplit(X,N,N), but the formula is simplified so
% it is much more efficient

d = double(bitand(x,2^n-1)) - double(bitget(x,n)).*2^n;
