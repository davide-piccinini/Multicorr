function vout = mcorr_signal_taper(vin,n)
%usage : vout=taper(vin,n)
%Applica un taper di dimensioni n campioni
%agli estremi del vettore vin.
%Davide Piccinini - Aprile 2000


[r,c]=size(vin);

if r < c
   vin=vin';
end

%t=taper (n,1);
t=hanning(n*2);
t=t;
ti=t(1:n);
tf=t(n+1:n*2);
l=length(vin);
vout=vin;
vin=vin;
vout(1:n)=vin(1:n).*ti;
vout(l-(n-1):l)=vin(l-(n-1):l).*tf;
