function c = mcorr_signal_filtrax(S,L,H,Fs);
%c = filtrax(S,L,H,Fs) Filtra un vettore S, usando un filtro di tipo Butterworth
%del 4 ordine.
%S=segnale in entrata
%H=limite superiore del filtro (Hz)
%L=limite inferiore del filtro (Hz)
%Fs=numero di campioni/secondo (cps)
%Davide Piccinini -Dicembre 1999-

if nargin < 4
    error('Usage: c=filtrax(S,L,H,Fs)');
end
if L & H ~=0
    [b,a]=butter(4,[L H]*2/Fs);
    %[H,w]=freqz(b,a,1024);
    c=filtfilt (b,a,S);
else
    c=S;
end
