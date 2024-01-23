function PICK=mcorr_io_get_picks(infile)

% PICK=mcorr_io_get_picks(infile);
% Read in sac file and buil a matlab structure with P and S arrivals

A=rsac(infile);
NZYEAR=lh(A,'NZYEAR');
NZJDAY=lh(A,'NZJDAY');
NZHOUR=lh(A,'NZHOUR');
NZMIN=lh(A,'NZMIN');
NZSEC=lh(A,'NZSEC');
NZMSEC=lh(A,'NZMSEC');
B     =lh(A,'B');
PT=lh(A,'A');
ST=lh(A,'T0');
[Y, D, M ]=jul2day(NZYEAR,NZJDAY);
%SACTIME=datenum(Y,M,D,NZHOUR,NZMIN,NZSEC+(NZMSEC/1000))+(B/86400);
SACTIME=datetime(Y,M,D,NZHOUR,NZMIN,NZSEC+(NZMSEC/1000))+(B/86400);
p=0;s=0;
if PT > 0
PARRIVAL=SACTIME+(PT/86400);
p=1;
end
if ST > 0
SARRIVAL=SACTIME+(ST/86400);
s=1;
end
PW=lh(A,'KA');PW=str2double(PW(4));
SW=lh(A,'KA');SW=str2double(SW(4));
if PW==0
Punc=0.01;
elseif PW==1
Punc=0.05;
elseif PW==2
Punc=0.1;
elseif PW==3
Punc=0.2;
end
if SW==0
Sunc=0.01;
elseif SW==1
Sunc=0.05;
elseif SW==2
Sunc=0.1;
elseif SW==3
Sunc=0.2;
end
PICK.sta=[];
PICK.net=[];
PICK.P=[];
PICK.Punc=[];
PICK.S=[];
PICK.Sunc=[];
PICK.sta=deblank(lh(A,'KSTNM'));
PICK.net=deblank(lh(A,'KNETWK'));
if p==1
PICK.P=PARRIVAL;
PICK.Punc=Punc;
end
if s==1
PICK.S=SARRIVAL;
PICK.Sunc=Sunc;
end 

PICK.data=A(:,2);
PICK.time=(A(:,1)/86400)+SACTIME;
PICK.DT  =lh(A,'DELTA');
PICK.CHN =deblank(lh(A,'KCMPNM'));




